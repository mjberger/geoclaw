# SGN Boussinesq — PETSc-free OpenMP linear solver (`bouss_solver = 1`)

This directory implements the implicit SGN (Serre–Green–Naghdi) Boussinesq update
for GeoClaw. Alongside the existing PETSc path it now has a **self-contained,
OpenMP-threaded linear solver** that needs neither PETSc nor MPI.

## Solver selection (`bouss.data` / `setrun.py`)

- `bouss_solver = 1` — OpenMP GMRES + directional block Gauss–Seidel preconditioner (this work)
- `bouss_solver = 2` — Pardiso (expired / unsupported)
- `bouss_solver = 3` — PETSc (retained; unchanged)

`isolver = bouss_solver` throughout the code.

## The linear system

The implicit update solves a sparse system stored in **CSR (compressed row)**
format with a **2×2 block** structure (`MatSetBlockSize 2` on the PETSc side).
Unknowns are interleaved per cell: `u0, v0, u1, v1, …` (even rows = x-momentum
`u`, odd rows = y-momentum `v`). Globally cells are numbered i-fastest.

```
A = [ I - D11    -D12   ]
    [  -D21    I - D22   ]
```

- **`I - D11`** (u-block) couples `U(i,j) ↔ U(i±1,j)` only ⇒ a direct sum of
  independent **1-D tridiagonals in x**.
- **`I - D22`** (v-block) couples `V(i,j) ↔ V(i,j±1)` only ⇒ independent
  **tridiagonals in y**.
- `-D12`, `-D21` are the u↔v cross-coupling blocks.

## How the solver works

1. `prepBuildSparseMatrixSGNcrs` assembles the CSR matrix for the *union* of
   Bouss grids at the level (threaded over grids), then `compressOut` compacts it.
2. `build_line_decomp` (called twice — u-lines parity 0, v-lines parity 1) walks
   the matrix graph and extracts the maximal tridiagonal **lines**. Each field
   dof lands in **exactly one line** (a true partition), which is what makes the
   threaded line solve race-free. The two decompositions of a level are
   **cached** in `line_cache(level)` together with the structure they were built
   from (rowPtr, cols, and which same-field off-diagonal entries were exactly
   zero, since those are not edges). `update_line_decomp` compares the current
   matrix with that structure in one threaded pass before every solve; if it is
   identical, only `diag/sub/sup` are refreshed through CSR positions saved at
   build time (`pdiag/psub/psup`), which is exactly what a rebuild would give.
   Any difference rebuilds both decompositions. The walk itself is still serial.
3. `gmres_solve` runs **right-preconditioned restarted GMRES(50)**, `rtol 1e-9`,
   up to 4 cycles, starting from the previous solve's psi (see "Nonzero initial
   guess"). The new Krylov vector is orthogonalized by **classical
   Gram-Schmidt in two threaded block passes** (`cgs_orth`: all j dot products in
   one sweep, all j updates in one sweep, re-orthogonalized only when the update
   removed more than a factor `cgs_reorth_tol = 1e-3` of the norm). Modified
   Gram-Schmidt is kept for comparison; see "Orthogonalization modes" below.
4. Preconditioner = **directional block Gauss–Seidel**: solve the u-block
   (tridiagonal lines in x, *exact* Thomas), update, solve the v-block. Forward
   sweep by default; symmetric (both cross blocks) available via
   `bgs_symmetric = .true.`.

## New source files

| File | Role |
|------|------|
| `bouss_tridiag_module.f90` | line decomposition (`build_line_decomp`), threaded Thomas line solve (`apply_line_solve`/`solve_one_path`), block-GS preconditioner (`apply_block_gs`), CSR sub-matvec |
| `bouss_gmres_module.f90` | OpenMP GMRES (`gmres_solve`), CSR matvec, threaded dot/axpy/scale |
| `test_tridiag.f90`, `test_gmres.f90`, `test_compress.f90` | standalone unit tests |

Modified: `implicit_update_bouss_2Calls.f90`, `prepBuildSparseMatrixSGNcrs.f90`,
`compressOut.f`, `amr2.f90`, `amr_module.f90`, `stst1.f`, `bouss_module.f90`,
`Makefile.bouss`.

## Building and running

- Builds **with or without `-fopenmp`** (serial fallback runs one thread).
- Builds **with or without `HAVE_PETSC`** — `bouss_solver = 1` needs neither PETSc
  nor MPI. Run single-process with `OMP_NUM_THREADS` set.
- PETSc-free `-O2` builds: `Makefile.nopetsc` in `examples/bouss/radial_flat`,
  `crater_westport` and `tohoku_sgn_hawaii` (only `CLAW` and `OMP_NUM_THREADS`
  needed): `make -f Makefile.nopetsc new; rm -f .data; make -f Makefile.nopetsc .data;
  make -f Makefile.nopetsc output OUTDIR=_output_xxx`. crater's `setrun.py` selects
  `bouss_solver = 3`; edit `bouss.data` to 1 after `make .data` (and `rm -f .data`
  first, or `make .data` is a no-op).
- Unit tests (not in the Makefile — compile by hand):
  ```
  gfortran -O2 -fopenmp bouss_tridiag_module.f90 bouss_gmres_module.f90 test_gmres.f90 -o t && ./t
  ```
  Tests are **bit-identical at 1 and N threads** — use them to check any change.

## Performance (crater_westport, 10 threads)

- `bouss_solver = 1`: **~219 s** linear solve — faster than PETSc redundant-LU
  (~230 s) and GAMG (~352 s), PETSc-free.
- Best with **`max1d = 60`**: more/smaller grids ⇒ better grid-parallelism in the
  (grid-threaded) matrix build. `compressOut` and `line-decomp` are per-level and
  grid-count independent.
- Diagnostics printed at end of run (`bouss_solver=1`): `line-decomp` and
  `compressOut` timers, the longest tridiagonal line seen (`maxLineLen`), how
  many line decompositions were rebuilt vs reused, the orthogonalization mode
  and re-orthogonalization count, the doubleCheck statistics, the wall time
  per Krylov phase (preconditioner / matvec / orthogonalization), the
  initial-guess statistics and the solves-by-kind table (see below).

## Orthogonalization modes (Nov 2026)

Selected at run time by environment variables, read once in `gmres_solve`:

- `BOUSS_GMRES_ORTH` = **3** (default) classical Gram-Schmidt, block passes;
  **1** modified Gram-Schmidt, one OpenMP region per dot/axpy (the original);
  **2** modified Gram-Schmidt serial (no OpenMP in the orthogonalization).
- `BOUSS_CGS_REORTH` = re-orthogonalization threshold for mode 3 (default `1e-3`;
  `0.7071` gives the DGKS / PETSc criterion).
- `doubleCheck` (module flag, default `.true.`; `BOUSS_GMRES_DOUBLECHECK=0` turns
  it off): a convergence claim from the Givens residual estimate is verified
  against the true residual ||b − A x|| before `gmres_solve` returns (one extra
  matvec per solve); if it fails, iteration continues with another restart
  cycle. The end-of-run summary prints the number of failed verifications and
  the worst |estimate − true| seen, so the stopping rule no longer depends on
  the orthogonalization keeping the basis orthogonal. Measured cost on
  radial_flat at `-O2`: none visible (68.4 s on vs 70.0 s off, inside the
  run-to-run noise); 0 failed verifications, max |est − true| = 1.6e-16.

radial_flat (200×200, 2 levels, 8 threads, M4 Max), linear-solve wall seconds;
iteration counts identical in every mode (11.00 / 19.07 on levels 1 / 2):

| build                        | MGS threaded | MGS serial | CGS 1e-3 | CGS DGKS |
|------------------------------|-------------:|-----------:|---------:|---------:|
| `-O -fbounds-check` (Makefile) | 111.7      | 235.5      | **91.6** | 119.1    |
| `-O2`, no bounds check        | 95.6        | 164.2      | **68.4** | 88.3     |

crater_westport (real bathymetry, 7 levels, 8 threads, `-O2`, `bouss_solver=1`):
MGS 427.1 s → CGS 368.4 s linear solve (−14%), total run 628.8 → 568.1 s (−10%);
iteration counts identical on all 7 levels (3.0 / 5.0 / 7.5 / 10.4 / 12.7 / 6.9 / 6.6);
224 re-orthogonalizations in ~143k iterations (0.16%); 0 failed doubleChecks.

tohoku_sgn_hawaii (real bathymetry, 6 levels, 8 threads, `-O2`, PETSc-free):
MGS 2743 s → CGS 2519 s linear solve (−8%; −10% excluding the line decompositions,
which were unchanged at ~700 s), total run 4144 → 3873 s (−7%); iteration counts
identical (2.0 / 3.0 / 4.7 / 9.0 / 4.2 / 8.9). The gain is smaller than on crater
because solves average only 2–9 iterations, so orthogonalization is a small share
of the solve. 15288 re-orthogonalizations in ~680k iterations (2.2%): they fire on
the *last* iteration of a solve, when the new vector is almost inside the Krylov
space (the same event as convergence), so the count scales with the number of
solves, not with any orthogonality problem; 0 failed doubleChecks, max
|est − true| = 1e-16.

Run-to-run noise between sessions is ~7% (the same MGS build measured 103.9 s in
October). A verification run (true residual recomputed after every solve) with
CGS 1e-3 gave max |est − true| = 3e-16 and max true relative residual 9.99e-10
over 2177 solves, so no re-orthogonalization was ever needed. The `-O2` build is what `xgeoclaw` in radial_flat currently is; a plain
`make new` there goes back to the Makefile's `-O -fbounds-check`.

## Line decomposition cache (Nov 2026)

`update_line_decomp` keeps both decompositions of a level in `line_cache(level)`
and only refreshes their values when the matrix structure is exactly the one they
were built from (see "How the solver works", item 2). Measured with CGS, 8 threads,
`-O2`, PETSc-free; "before" is the same-day CGS run of each case:

| case (levels)         | solves | rebuilt | reused | line-decomp s   | linear solve s           |
|-----------------------|-------:|--------:|-------:|----------------:|-------------------------:|
| radial_flat (2)       |   2179 |     172 |   2007 | 11.9 → **2.7**  | 68.4 → **57.5** (−16%)   |
| crater_westport (7)   |  15488 |     482 |  15006 | 114.1 → **15.0**| 368.4 → **250.5** (−32%) |
| tohoku_sgn_hawaii (6) |  88098 |     779 |  87319 | 699.0 → **91.8**| 2518.6 → **1644.1** (−35%)|

Iteration counts unchanged in all three; 0 failed doubleChecks. A rebuild happens
only after regridding touched that level (1 per 13 / 32 / 113 solves). The honest
cache gain is the line-decomp drop, 13% / 27% / 24% of the old solve time; the rest
of the measured change tracks the solver-independent timers (Stepgrid was 2% / 9% /
20% faster in the "after" runs, i.e. the earlier crater and tohoku runs were made
while the laptop was also in use). **Always compare a solver-independent timer
before attributing a gain.** What remains of `line-decomp` is the threaded
structure comparison + value refresh (~1 ms per solve on tohoku) plus the rebuilds,
whose walk is still serial. Outputs: `_output_O2_cache` in each example directory.

## Nonzero initial guess (Nov 2026)

Every solve starts from the psi of the previous solve on that level instead of
from zero (module flag `gmres_use_guess`, default on; `BOUSS_GMRES_GUESS=0`
turns it off).
`buildSparseMatrixSGNcrs` loads it into `soln` from components 4,5 of the
old-time storage (`qold`, which is `q` on the finest level) while it numbers
the cells, so regridding needs no special handling (new grids carry
interpolated psi; cells that revert to SWE get 0). `gmres_solve` measures
||b − A x0||/||b|| with the first residual evaluation it makes anyway and, if
the guess is worse than zero, falls back to x0 = 0. That guard is needed: on
nearly quiescent grids ||b|| is at round-off level while the stored psi is
not (858 of 15488 guesses rejected on crater, ratios up to 1e12). The stopping
test stays ||b − A x|| ≤ rtol ||b||, so accuracy is unchanged: 0 failed
doubleChecks, same max true residual, radial_flat final frame identical in
grid structure with values differing by 2e-10 relative.

| case (8 threads, `-O2`) | linear solve s, off → on | avg iterations per level, off → on |
|---|---|---|
| radial_flat (2 levels)     | 55.1 → **47.4** (−14%)  | 11.0 / 19.1 → 9.7 / 16.7 |
| crater_westport (7 levels) | 247.4 → **210.3** (−15%) | 3.0 / 5.0 / 7.5 / 10.4 / 12.7 / 6.9 / 6.6 → 3.0 / 4.6 / 7.1 / 8.8 / 9.6 / 6.6 / 6.0 |

The solves split into three kinds (`gmres_guess_kind`, set in `implicit_update`
from `doUpdate` and `levelBouss < mxnest`); radial_flat, guess on:

| kind | the guess is | solves | iters off → on | mean eff. ratio |
|---|---|---:|---:|---:|
| 1 coarse-level 1st solve | previous step's 2nd solve, same time level, but the fine-grid overlap was replaced by fine averages since | 341 | 11.00 → 9.36 | 0.033 |
| 2 coarse-level 2nd solve | this step's 1st solve, one SWE step stale | 341 | 11.00 → 10.00 | 0.23 |
| 3 finest level | previous step's solve, one step stale | 1497 | 19.10 → 16.65 | 0.064 |

A 10–30× smaller initial residual saves only 1.5–2.5 iterations because the
preconditioned GMRES contracts by ~0.15–0.35 per iteration, so this is close
to the ceiling for any initial guess. Re-orthogonalizations rise (224 → 1016
on crater) because more solves end on the nearly-invariant-Krylov-space event;
harmless. An earlier attempt at an initial guess (PETSc era) showed no gain;
a stopping test relative to the initial residual instead of ||b|| would
explain that.

Krylov phase timers, guess off: radial_flat 55 s = preconditioner 22.4 +
orthogonalization 19.9 + matvec 7.5 + rest 5 (line-decomp 2.6, vector ops);
crater 247 s = 113.6 + 66.6 + 37.0 + rest 30 (line-decomp 15). The block
Gauss–Seidel apply is the next target (fused vector passes in `apply_block_gs`
/ `sub_matvec`, one reciprocal per line in the Thomas solve).

## Hard-won lessons / gotchas (read before editing)

0. **DGKS fires on every GMRES iteration.** The textbook "re-orthogonalize if
   ||w|| dropped by more than 1/√2" test (PETSc's default too) assumes nothing
   about w. In GMRES on a well-preconditioned system the new vector A·M⁻¹v lies
   mostly in the span of the old ones, so the projection always removes most of
   its norm and CGS silently becomes CGS2 (36008 refinements in 36012 iterations
   here), slower than MGS. GMRES only needs orthogonality loss ≈ ε/threshold to
   stay far below `rtol`, hence the 1e-3 default. **Serial MGS is 2× slower**
   than threaded MGS: the dot/axpy regions are bandwidth-bound, not sync-bound,
   so "fewer forks" is not where the time goes; the CGS gain is cache reuse of
   the thread's chunk of w across all j columns plus 2 regions instead of 2j.

1. **OpenMP without a hard dependency.** Do **not** `use omp_lib`. Declare
   `omp_get_max_threads`/`omp_get_thread_num` as `integer` externals and call them
   **only on `!$` sentinel lines**, defaulting thread count/id to 1/0. (Matches
   `amr2`/`advanc`/`flagger`.) This is what lets the code build without `-fopenmp`.

2. **Persistent `save` work buffers.** The big arrays (GMRES Krylov basis `V`
   ≈ 80 MB, block-GS `Ay`/`rw`, Thomas `thom_cp`, `compressOut` scratch) are
   module-`save`, grown-as-needed, reused across solves. `gmres_solve` /
   `apply_block_gs` are entered **only on the master thread** (`implicit_update`
   is called serially per level), so `save` is thread-safe. Eliminating per-call
   allocation of `V` saved ~27 s on crater.

3. **Worker-thread stacks.** OpenMP *worker* threads get a small stack
   (`OMP_STACKSIZE`); only the *primary* thread has the big main stack. Automatic
   (stack) arrays sized by problem size overflow worker stacks on large cases even
   when a single-thread run is fine. That is why the Thomas scratch is on the
   **heap** (`thom_cp(maxline, nthreads)`), not an automatic array. Tridiagonal
   lines cross grid boundaries, so line length is **not** bounded by `max1d`.
   On **macOS**, `limit stacksize unlimited` is silently capped (~64 MB main
   thread), and `OMP_STACKSIZE` only affects workers — so bumping the stack may
   not "fix" an overflow.

4. **The line decomposition must stay a true partition** (each dof in exactly one
   line). The parallel line solve writes disjoint dofs and reads read-only shared
   data; correctness depends on the partition. `build_line_decomp` guarantees it
   via a `visited` flag (max degree 2 asserted).

5. **Debugging NaNs.** A *quiet* NaN from uninitialized memory propagates silently
   through `+ − × ÷` and only traps at an **ordered comparison** under
   `-ffpe-trap=invalid` — so the trap fires where it's *used*, not created. Use
   **`-finit-real=snan`** (with `-ffpe-trap=invalid`) to trap at first *touch* and
   find the source. (A global-case NaN was tracked this way — it was in the
   coarse ghost-cell fill, **not** the solver.)

6. **Debugging on macOS.** `gfortran -fbacktrace` prints only raw addresses, and
   ASLR makes them unsymbolizable after the fact. Use **lldb**
   (`lldb --batch -o run -k bt -- ./xgeoclaw`) for symbolized traces.

## Provenance

GMRES (Saad & Schultz 1986), the Thomas tridiagonal solve (1949), graph traversal
for line extraction, and stream compaction are all standard algorithms, written
from scratch here — not copied from PETSc, SPARSKIT, Numerical Recipes, etc.
