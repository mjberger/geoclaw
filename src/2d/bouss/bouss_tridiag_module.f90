! =====================================================================
! bouss_tridiag_module
!
! Tridiagonal line preconditioner for the SGN Boussinesq linear system,
! used by the OpenMP GMRES in bouss_gmres_module (bouss_solver = 1).
! The 2x2-block matrix (interleaved u0,v0,u1,v1,...) has
! diagonal blocks that are each a direct sum of INDEPENDENT 1D tridiagonal
! systems:
!    field 0 (u, even 0-based rows): I-D11 couples only U(i,j)<->U(i+-1,j)
!    field 1 (v, odd  0-based rows): I-D22 couples only V(i,j)<->V(i,j+-1)
!
! This module (a) extracts those tridiagonal lines directly from the
! assembled compressed CRS matrix -- graph-based, so it is agnostic to
! grid geometry, patch boundaries, BC folding and revert-to-SWE cells --
! and (b) solves them with the Thomas algorithm, threaded over lines with
! OpenMP (the lines are independent, so this needs no communication).
!
! Field membership is by parity of the 0-based row index: even = u, odd = v
! (consistent with MatSetBlockSize(J,2), field 0 = component 0).
!
! Standalone: depends only on OpenMP, no PETSc / amr_module / bouss_module,
! so it can be unit-tested with plain gfortran before integration.
! =====================================================================
module bouss_tridiag_module

    implicit none
    private

    public :: line_decomp_t
    public :: build_line_decomp
    public :: apply_line_solve
    public :: free_line_decomp
    public :: sub_matvec
    public :: apply_block_gs
    public :: validate_line_solve
    public :: line_cache_t, line_cache, update_line_decomp, ensure_line_cache
    public :: refresh_line_values

    ! set .false. to disable the (once-per-level) round-trip validation check
    logical, public :: tridiag_validate = .true.

    ! .false. (default) = FORWARD block Gauss-Seidel (solve u, then v using A10;
    !   ignores A01).  Matches PETSc multiplicative fieldsplit and is cheaper
    !   (2 block solves + 1 cross matvec per apply).  .true. = SYMMETRIC sweep
    !   (also re-solves u using A01 y_v): stronger for very strong u<->v coupling
    !   but ~1.5x the apply cost.  Forward converges the tested cases (crater,
    !   radial); flip to .true. only if a case stalls with forward.
    logical, public :: bgs_symmetric = .false.

    ! persistent work buffers for apply_block_gs, grown as needed and reused
    ! across applies (avoid re-allocating every GMRES iteration).  apply_block_gs
    ! is entered on the master thread only, so save is thread-safe.
    real(kind=8), allocatable, save :: Ay(:), rw(:)
    integer, save :: bgs_ntot = -1

    ! Decomposition of one field's diagonal block into tridiagonal lines.
    ! Storage is CSR-over-paths: path p occupies ordered slots
    ! pathStart(p) .. pathStart(p+1)-1 in the dof/sub/diag/sup arrays.
    type :: line_decomp_t
        integer :: ntot   = 0          ! total dofs in the full 2N system
        integer :: parity = -1         ! 0 = u (even rows), 1 = v (odd rows)
        integer :: nfield = 0          ! number of dofs in this field (= N)
        integer :: npaths = 0          ! number of tridiagonal lines
        integer, allocatable :: pathStart(:)  ! (npaths+1)
        integer, allocatable :: dof(:)        ! (nfield) 0-based row index, path-ordered
        real(kind=8), allocatable :: sub(:)   ! (nfield) coupling to previous dof in path
        real(kind=8), allocatable :: diag(:)  ! (nfield) diagonal
        real(kind=8), allocatable :: sup(:)   ! (nfield) coupling to next dof in path
        ! CSR positions the three coefficient arrays were read from (-1 = not
        ! stored, value 0), so they can be refreshed without rebuilding paths
        integer, allocatable :: pdiag(:), psub(:), psup(:)   ! (nfield)
        integer :: maxline = 0                ! longest path (for scratch sizing)
        logical :: has_cycle = .false.        ! true if a periodic (cyclic) line was found
    end type line_decomp_t

    ! Per-level cache of the u/v decompositions together with the matrix
    ! structure they were built from (rowPtr, cols, and which same-field
    ! off-diagonal entries were exactly zero).  See update_line_decomp.
    type :: line_cache_t
        logical :: valid = .false.
        integer :: ntot = -1, nnz = -1
        integer, allocatable :: rowPtr(:)     ! (0:ntot)
        integer, allocatable :: cols(:)       ! (0:nnz-1)
        logical, allocatable :: zmask(:)      ! (0:nnz-1) structural-zero flags
        type(line_decomp_t) :: D0, D1
    end type line_cache_t
    type(line_cache_t), allocatable, save :: line_cache(:)    ! indexed by level
    integer(kind=8), public, save :: linecache_nreuse = 0, linecache_nrebuild = 0

    ! persistent per-thread Thomas scratch (maxline x nthreads), grown as needed,
    ! reused across applies.  Keeps the O(line-length) work off the WORKER-thread
    ! stack (automatic arrays there overflow the small OMP stack on big grids;
    ! that was a real crash).  apply_line_solve is called serially (master
    ! thread) so growing it is race-free; the parallel loop only reads its shape
    ! and writes disjoint columns.
    real(kind=8), allocatable, save :: thom_cp(:,:)

contains

    ! -----------------------------------------------------------------
    ! Look up A(r,c) from compressed CRS (0-based row r, column c).
    ! Returns 0 if not stored.
    ! -----------------------------------------------------------------
    pure function get_entry(rowPtr, cols, vals, nnz, ntot, r, c) result(v)
        integer,      intent(in) :: nnz, ntot, r, c
        integer,      intent(in) :: rowPtr(0:ntot), cols(0:nnz-1)
        real(kind=8), intent(in) :: vals(0:nnz-1)
        real(kind=8) :: v
        integer :: p
        v = 0.d0
        do p = rowPtr(r), rowPtr(r+1)-1
            if (cols(p) == c) then
                v = vals(p)
                return
            endif
        end do
    end function get_entry

    ! index (1 or 2) of the same-parity edge from a dof to targetf (field-local)
    ! among the (<=2) off-diagonal entries captured from that dof's row, else 0.
    pure function edge_index(edof, ecnt, targetf) result(k)
        integer, intent(in) :: edof(2), ecnt, targetf
        integer :: k, kk
        k = 0
        do kk = 1, ecnt
            if (edof(kk) == targetf) then
                k = kk
                return
            endif
        end do
    end function edge_index

    ! -----------------------------------------------------------------
    ! Build the tridiagonal-line decomposition of one field's diagonal
    ! block from the assembled compressed CRS of the full 2N system.
    !
    !   rowPtr(0:ntot), cols(0:nnz-1), vals(0:nnz-1)  : 0-based CRS
    !   ntot   : number of rows/cols of the full system (= 2N)
    !   parity : 0 for u (even rows), 1 for v (odd rows)
    !   D      : output decomposition
    !
    ! Method: within the chosen field, connect two dofs if EITHER
    ! A(r,c) or A(c,r) is a stored nonzero (symmetric closure) -- this
    ! keeps revert cells (whose own row lost its off-diagonals) attached
    ! to their neighbours.  For a genuine 1D tridiagonal operator every
    ! node then has degree <= 2, so the field is a union of simple paths
    ! (open lines) plus isolated nodes.  Cyclic lines (periodic BCs) are
    ! detected and flagged (not yet solved specially).
    ! -----------------------------------------------------------------
    subroutine build_line_decomp(rowPtr, cols, vals, nnz, ntot, parity, D)
        integer,      intent(in) :: nnz, ntot, parity
        integer,      intent(in) :: rowPtr(0:ntot), cols(0:nnz-1)
        real(kind=8), intent(in) :: vals(0:nnz-1)
        type(line_decomp_t), intent(out) :: D

        integer :: nfield, f, r, c, p, cf
        integer, allocatable :: dof2f(:)     ! full dof -> field-local (1-based), else -1
        integer, allocatable :: f2dof(:)     ! field-local -> full dof (0-based)
        integer, allocatable :: nbr(:,:)     ! (2,nfield) field-local neighbours, -1 if none
        integer, allocatable :: deg(:)       ! (nfield) degree
        logical, allocatable :: visited(:)
        integer :: np, slot, cur, prev, nxt, startf, prevf, nextf, kk
        ! values captured during the adjacency scan (avoid a 2nd row scan):
        real(kind=8), allocatable :: diagval(:)     ! (nfield) A(r,r)
        integer,      allocatable :: diagpos(:)     ! (nfield) CSR position of A(r,r), -1 if none
        integer,      allocatable :: edge_dof(:,:)  ! (2,nfield) same-parity off-diag nbrs (field-local)
        real(kind=8), allocatable :: edge_val(:,:)  ! (2,nfield) A(r, edge_dof) from row r
        integer,      allocatable :: edge_pos(:,:)  ! (2,nfield) CSR position of edge_val
        integer,      allocatable :: edge_cnt(:)    ! (nfield) # off-diag entries in row r (0,1,2)

        ! ---- map field dofs <-> field-local indices ----
        nfield = 0
        do r = parity, ntot-1, 2
            nfield = nfield + 1
        end do

        allocate(dof2f(0:ntot-1)); dof2f = -1
        allocate(f2dof(nfield))
        f = 0
        do r = parity, ntot-1, 2
            f = f + 1
            f2dof(f) = r
            dof2f(r) = f
        end do

        ! ---- one pass per field row: build symmetric-closure adjacency AND
        ! capture the diagonal + directional off-diagonal values, so the
        ! coefficient-extraction pass below needs no second row scan. ----
        allocate(nbr(2,nfield)); nbr = -1
        allocate(deg(nfield));   deg = 0
        allocate(diagval(nfield));    diagval  = 0.d0
        allocate(diagpos(nfield));    diagpos  = -1
        allocate(edge_dof(2,nfield)); edge_dof = -1
        allocate(edge_val(2,nfield)); edge_val = 0.d0
        allocate(edge_pos(2,nfield)); edge_pos = -1
        allocate(edge_cnt(nfield));   edge_cnt = 0

        do f = 1, nfield
            r = f2dof(f)
            do p = rowPtr(r), rowPtr(r+1)-1
                c = cols(p)
                if (c < 0 .or. c > ntot-1) cycle
                if (mod(c,2) /= parity) cycle      ! other field (cross-coupling), skip
                if (c == r) then
                    diagval(f) = vals(p)           ! diagonal A(r,r)
                    diagpos(f) = p
                    cycle
                endif
                if (vals(p) == 0.d0) cycle         ! explicit zero, not an edge
                cf = dof2f(c)
                call add_nbr(nbr, deg, nfield, f, cf)  ! symmetric topological edge
                if (edge_cnt(f) < 2) then          ! directional value A(r,c)
                    edge_cnt(f) = edge_cnt(f) + 1
                    edge_dof(edge_cnt(f), f) = cf
                    edge_val(edge_cnt(f), f) = vals(p)
                    edge_pos(edge_cnt(f), f) = p
                endif
            end do
        end do

        ! ---- walk paths ----
        allocate(D%pathStart(nfield+1))   ! at most nfield paths (all singletons)
        allocate(D%dof(nfield))
        allocate(visited(nfield)); visited = .false.

        np   = 0
        slot = 0
        D%has_cycle = .false.

        ! First: every path that has an endpoint (degree 0 or 1).
        do startf = 1, nfield
            if (visited(startf)) cycle
            if (deg(startf) > 1) cycle          ! interior node, handle in walk from an end
            np = np + 1
            D%pathStart(np) = slot + 1
            prev = -1
            cur  = startf
            do
                visited(cur) = .true.
                slot = slot + 1
                D%dof(slot) = f2dof(cur)
                nxt = other_nbr(nbr, cur, prev)
                if (nxt == -1) exit
                if (visited(nxt)) exit           ! safety
                prev = cur
                cur  = nxt
            end do
        end do

        ! Any node still unvisited is part of a cycle (all degree 2).
        do startf = 1, nfield
            if (visited(startf)) cycle
            D%has_cycle = .true.
            np = np + 1
            D%pathStart(np) = slot + 1
            prev = -1
            cur  = startf
            do
                if (visited(cur)) exit
                visited(cur) = .true.
                slot = slot + 1
                D%dof(slot) = f2dof(cur)
                nxt = other_nbr(nbr, cur, prev)
                if (nxt == -1) exit
                prev = cur
                cur  = nxt
            end do
        end do

        D%pathStart(np+1) = slot + 1
        D%npaths = np
        D%nfield = nfield
        D%ntot   = ntot
        D%parity = parity
        D%maxline = 0
        do np = 1, D%npaths
            D%maxline = max(D%maxline, D%pathStart(np+1) - D%pathStart(np))
        end do

        ! ---- extract tridiagonal coefficients along each ordered path ----
        ! uses the values captured in the adjacency scan (no 2nd row scan)
        allocate(D%sub(nfield), D%diag(nfield), D%sup(nfield))
        allocate(D%psub(nfield), D%pdiag(nfield), D%psup(nfield))
        D%sub = 0.d0; D%diag = 0.d0; D%sup = 0.d0
        D%psub = -1;  D%pdiag = -1;  D%psup = -1
        do np = 1, D%npaths
            do slot = D%pathStart(np), D%pathStart(np+1)-1
                f = dof2f(D%dof(slot))
                D%diag(slot)  = diagval(f)
                D%pdiag(slot) = diagpos(f)
                if (slot > D%pathStart(np)) then
                    prevf = dof2f(D%dof(slot-1))
                    kk = edge_index(edge_dof(:,f), edge_cnt(f), prevf)
                    if (kk > 0) then
                        D%sub(slot) = edge_val(kk,f);  D%psub(slot) = edge_pos(kk,f)
                    endif
                endif
                if (slot < D%pathStart(np+1)-1) then
                    nextf = dof2f(D%dof(slot+1))
                    kk = edge_index(edge_dof(:,f), edge_cnt(f), nextf)
                    if (kk > 0) then
                        D%sup(slot) = edge_val(kk,f);  D%psup(slot) = edge_pos(kk,f)
                    endif
                endif
            end do
        end do

        deallocate(dof2f, f2dof, nbr, deg, visited)
        deallocate(diagval, diagpos, edge_dof, edge_val, edge_pos, edge_cnt)
    end subroutine build_line_decomp

    ! insert field-local neighbour j into node i's list (and i into j's),
    ! avoiding duplicates; abort if a node would exceed degree 2.
    subroutine add_nbr(nbr, deg, nfield, i, j)
        integer, intent(in)    :: nfield, i, j
        integer, intent(inout) :: nbr(2,nfield), deg(nfield)
        call add_one(nbr, deg, nfield, i, j)
        call add_one(nbr, deg, nfield, j, i)
    end subroutine add_nbr

    subroutine add_one(nbr, deg, nfield, i, j)
        integer, intent(in)    :: nfield, i, j
        integer, intent(inout) :: nbr(2,nfield), deg(nfield)
        if (nbr(1,i) == j .or. nbr(2,i) == j) return   ! already present
        if (deg(i) == 0) then
            nbr(1,i) = j; deg(i) = 1
        else if (deg(i) == 1) then
            nbr(2,i) = j; deg(i) = 2
        else
            write(*,*) 'bouss_tridiag: ERROR node ',i, &
                       ' has degree > 2 (block is not tridiagonal)'
            stop 1
        endif
    end subroutine add_one

    ! given node cur and the node we came from (prev, -1 if none),
    ! return the other neighbour, or -1 if none.
    pure function other_nbr(nbr, cur, prev) result(nxt)
        integer, intent(in) :: nbr(:,:), cur, prev
        integer :: nxt
        if (nbr(1,cur) /= prev .and. nbr(1,cur) /= -1) then
            nxt = nbr(1,cur)
        else if (nbr(2,cur) /= prev .and. nbr(2,cur) /= -1) then
            nxt = nbr(2,cur)
        else
            nxt = -1
        endif
    end function other_nbr

    ! -----------------------------------------------------------------
    ! Solve  A_block * x = rhs  for one field, exactly, by Thomas on each
    ! independent line, threaded over lines.  rhs and x are full-length
    ! (0:ntot-1); only the field's dofs of x are written.
    ! -----------------------------------------------------------------
    ! x is intent(inout): only this field's dofs are written, so entries
    ! belonging to the other field are left untouched (relied on by
    ! apply_block_gs, which solves u then v into the same vector).
    subroutine apply_line_solve(D, rhs, x)
        type(line_decomp_t), intent(in)    :: D
        real(kind=8),        intent(in)    :: rhs(0:D%ntot-1)
        real(kind=8),        intent(inout) :: x(0:D%ntot-1)
        integer :: p, tid, nth
        integer :: omp_get_max_threads, omp_get_thread_num   ! external (only
        !          called on !$ lines; harmless declaration without OpenMP)

        ! grow the persistent per-thread scratch if needed (serial context here).
        ! Without OpenMP nth defaults to 1 (single column of scratch).
        nth = 1
        !$ nth = omp_get_max_threads()
        if (.not. allocated(thom_cp)) then
            allocate(thom_cp(max(D%maxline,1), nth))
        else if (size(thom_cp,1) < D%maxline .or. size(thom_cp,2) < nth) then
            deallocate(thom_cp)
            allocate(thom_cp(max(D%maxline,1), nth))
        endif

        !$omp parallel do schedule(dynamic,16) default(shared) private(p,tid)
        do p = 1, D%npaths
            tid = 1
            !$ tid = omp_get_thread_num() + 1
            call solve_one_path(D, p, rhs, x, thom_cp(:,tid))
        end do
        !$omp end parallel do
    end subroutine apply_line_solve

    ! Thomas solve of one tridiagonal line, reading coefficients directly from
    ! D and rhs (no per-call copies) and writing into x at the line's dofs.
    ! cp is thread-private scratch (length >= line length); the forward pass
    ! stores dp in x, the back pass overwrites it in place.  No large automatic
    ! (stack) arrays -> safe on small OpenMP worker stacks.
    subroutine solve_one_path(D, p, rhs, x, cp)
        type(line_decomp_t), intent(in)    :: D
        integer,             intent(in)    :: p
        real(kind=8),        intent(in)    :: rhs(0:D%ntot-1)
        real(kind=8),        intent(inout) :: x(0:D%ntot-1)
        real(kind=8),        intent(inout) :: cp(:)
        integer :: off, n, m
        real(kind=8) :: denom, ai, bi

        off = D%pathStart(p) - 1
        n   = D%pathStart(p+1) - D%pathStart(p)

        ! forward elimination (store dp in x at the line's dof positions)
        bi = D%diag(off+1)
        cp(1) = D%sup(off+1) / bi
        x(D%dof(off+1)) = rhs(D%dof(off+1)) / bi
        do m = 2, n
            ai    = D%sub(off+m)
            denom = D%diag(off+m) - ai*cp(m-1)
            cp(m) = D%sup(off+m) / denom
            x(D%dof(off+m)) = (rhs(D%dof(off+m)) - ai*x(D%dof(off+m-1))) / denom
        end do

        ! back substitution (in place, same thread, disjoint from other paths)
        do m = n-1, 1, -1
            x(D%dof(off+m)) = x(D%dof(off+m)) - cp(m)*x(D%dof(off+m+1))
        end do
    end subroutine solve_one_path

    subroutine free_line_decomp(D)
        type(line_decomp_t), intent(inout) :: D
        if (allocated(D%pathStart)) deallocate(D%pathStart)
        if (allocated(D%dof))       deallocate(D%dof)
        if (allocated(D%sub))       deallocate(D%sub)
        if (allocated(D%diag))      deallocate(D%diag)
        if (allocated(D%sup))       deallocate(D%sup)
        if (allocated(D%pdiag))     deallocate(D%pdiag)
        if (allocated(D%psub))      deallocate(D%psub)
        if (allocated(D%psup))      deallocate(D%psup)
        D%ntot = 0; D%parity = -1; D%nfield = 0; D%npaths = 0
        D%has_cycle = .false.
    end subroutine free_line_decomp

    ! -----------------------------------------------------------------
    ! Decomposition cache.
    !
    ! The u/v line decompositions of a level depend only on the matrix
    ! STRUCTURE: rowPtr, cols, and which same-field off-diagonal entries are
    ! exactly zero (build_line_decomp does not treat those as edges).  Between
    ! regrids the structure of a level is unchanged for many solves while the
    ! values change every solve.  update_line_decomp compares the current
    ! matrix against the structure the cached decompositions were built from,
    ! in one threaded pass over the index arrays; if identical, only the
    ! coefficient arrays are refreshed through the CSR positions saved at
    ! build time (also threaded), which gives exactly what a rebuild would.
    ! Any difference at all triggers the full rebuild.
    ! -----------------------------------------------------------------

    ! make sure line_cache(1:nlev) exists (entries are kept across growth)
    subroutine ensure_line_cache(nlev)
        integer, intent(in) :: nlev
        type(line_cache_t), allocatable :: tmp(:)
        integer :: n
        if (.not. allocated(line_cache)) then
            allocate(line_cache(max(nlev,1)))
        else if (size(line_cache) < nlev) then
            n = size(line_cache)
            allocate(tmp(nlev))
            tmp(1:n) = line_cache(1:n)
            call move_alloc(tmp, line_cache)
        endif
    end subroutine ensure_line_cache

    ! structural-zero flag of CSR entry (row r, column c, value v):
    ! same field, off-diagonal, exactly zero -> not an edge for build_line_decomp
    pure logical function szero(r, c, v)
        integer,      intent(in) :: r, c
        real(kind=8), intent(in) :: v
        szero = (c /= r) .and. (mod(c,2) == mod(r,2)) .and. (v == 0.d0)
    end function szero

    ! does the current matrix have exactly the structure cached in c?
    logical function same_structure(c, rowPtr, cols, vals, nnz, ntot) result(same)
        type(line_cache_t), intent(in) :: c
        integer,      intent(in) :: nnz, ntot, rowPtr(0:ntot), cols(0:nnz-1)
        real(kind=8), intent(in) :: vals(0:nnz-1)
        integer :: r, p
        logical :: ok
        same = .false.
        if (.not. c%valid) return
        if (c%ntot /= ntot .or. c%nnz /= nnz) return
        if (rowPtr(ntot) /= c%rowPtr(ntot)) return
        ok = .true.
        !$omp parallel do schedule(static) default(shared) private(r,p) reduction(.and.:ok)
        do r = 0, ntot-1
            if (rowPtr(r) /= c%rowPtr(r)) ok = .false.
            do p = rowPtr(r), rowPtr(r+1)-1
                if (cols(p) /= c%cols(p)) ok = .false.
                if (szero(r, cols(p), vals(p)) .neqv. c%zmask(p)) ok = .false.
            end do
        end do
        !$omp end parallel do
        same = ok
    end function same_structure

    ! remember the structure the decompositions in c were just built from
    subroutine store_structure(c, rowPtr, cols, vals, nnz, ntot)
        type(line_cache_t), intent(inout) :: c
        integer,      intent(in) :: nnz, ntot, rowPtr(0:ntot), cols(0:nnz-1)
        real(kind=8), intent(in) :: vals(0:nnz-1)
        integer :: r, p
        if (allocated(c%rowPtr)) deallocate(c%rowPtr, c%cols, c%zmask)
        allocate(c%rowPtr(0:ntot), c%cols(0:max(nnz,1)-1), c%zmask(0:max(nnz,1)-1))
        c%ntot = ntot; c%nnz = nnz
        c%rowPtr(ntot) = rowPtr(ntot)
        !$omp parallel do schedule(static) default(shared) private(r,p)
        do r = 0, ntot-1
            c%rowPtr(r) = rowPtr(r)
            do p = rowPtr(r), rowPtr(r+1)-1
                c%cols(p)  = cols(p)
                c%zmask(p) = szero(r, cols(p), vals(p))
            end do
        end do
        !$omp end parallel do
        c%valid = .true.
    end subroutine store_structure

    ! re-read diag/sub/sup from the current values through the saved positions
    subroutine refresh_line_values(D, vals, nnz)
        type(line_decomp_t), intent(inout) :: D
        integer,      intent(in) :: nnz
        real(kind=8), intent(in) :: vals(0:nnz-1)
        integer :: s
        !$omp parallel do schedule(static) default(shared) private(s)
        do s = 1, D%nfield
            if (D%pdiag(s) >= 0) then
                D%diag(s) = vals(D%pdiag(s))
            else
                D%diag(s) = 0.d0
            endif
            if (D%psub(s) >= 0) then
                D%sub(s) = vals(D%psub(s))
            else
                D%sub(s) = 0.d0
            endif
            if (D%psup(s) >= 0) then
                D%sup(s) = vals(D%psup(s))
            else
                D%sup(s) = 0.d0
            endif
        end do
        !$omp end parallel do
    end subroutine refresh_line_values

    ! Bring line_cache(level)%D0 / %D1 up to date with the current matrix:
    ! refresh the coefficients if the structure is unchanged, else rebuild.
    subroutine update_line_decomp(level, rowPtr, cols, vals, nnz, ntot, rebuilt)
        integer,      intent(in)  :: level, nnz, ntot, rowPtr(0:ntot), cols(0:nnz-1)
        real(kind=8), intent(in)  :: vals(0:nnz-1)
        logical,      intent(out) :: rebuilt
        call ensure_line_cache(level)
        if (same_structure(line_cache(level), rowPtr, cols, vals, nnz, ntot)) then
            call refresh_line_values(line_cache(level)%D0, vals, nnz)
            call refresh_line_values(line_cache(level)%D1, vals, nnz)
            rebuilt = .false.
            linecache_nreuse = linecache_nreuse + 1
        else
            call build_line_decomp(rowPtr, cols, vals, nnz, ntot, 0, line_cache(level)%D0)
            call build_line_decomp(rowPtr, cols, vals, nnz, ntot, 1, line_cache(level)%D1)
            call store_structure(line_cache(level), rowPtr, cols, vals, nnz, ntot)
            rebuilt = .true.
            linecache_nrebuild = linecache_nrebuild + 1
        endif
    end subroutine update_line_decomp

    ! -----------------------------------------------------------------
    ! Sub-block matvec:  b(r) = sum_{c: parity(c)=col_parity} A(r,c) x(c)
    ! over rows r with parity(r)=row_parity; b = 0 elsewhere.
    !   row_parity==col_parity : a diagonal block (A00 or A11)
    !   row_parity/=col_parity : a cross block   (A10 or A01)
    ! -----------------------------------------------------------------
    subroutine sub_matvec(rowPtr, cols, vals, nnz, ntot, row_parity, col_parity, x, b)
        integer,      intent(in)  :: nnz, ntot, row_parity, col_parity
        integer,      intent(in)  :: rowPtr(0:ntot), cols(0:nnz-1)
        real(kind=8), intent(in)  :: vals(0:nnz-1), x(0:ntot-1)
        real(kind=8), intent(out) :: b(0:ntot-1)
        integer :: r, p, c
        real(kind=8) :: s
        b = 0.d0
        !$omp parallel do schedule(static) default(shared) private(r,p,c,s)
        do r = row_parity, ntot-1, 2
            s = 0.d0
            do p = rowPtr(r), rowPtr(r+1)-1
                c = cols(p)
                if (mod(c,2) == col_parity) s = s + vals(p)*x(c)
            end do
            b(r) = s
        end do
        !$omp end parallel do
    end subroutine sub_matvec

    ! -----------------------------------------------------------------
    ! SYMMETRIC block Gauss-Seidel preconditioner apply (matrix-free):
    !     y_u = A00^{-1} r_u                     (forward: u)
    !     y_v = A11^{-1} ( r_v - A10 y_u )       (forward: v, uses A10)
    !     y_u = A00^{-1} ( r_u - A01 y_v )       (backward: re-solve u, uses A01)
    ! Unlike a forward-only sweep (which ignores A01), this uses BOTH cross
    ! blocks, so it stays effective when the u<->v coupling is strong (real
    ! bathymetry: the topographic -D12/-D21 terms).  When A01 ~ 0 it reduces
    ! to the forward sweep, so weak-coupling cases are unchanged.
    ! D0/D1 are the u/v line decompositions of the SAME matrix (rowPtr,cols,
    ! vals); block solves are the OpenMP Thomas line solves.
    ! -----------------------------------------------------------------
    subroutine apply_block_gs(D0, D1, rowPtr, cols, vals, nnz, ntot, r, y)
        type(line_decomp_t), intent(in)  :: D0, D1
        integer,      intent(in)  :: nnz, ntot
        integer,      intent(in)  :: rowPtr(0:ntot), cols(0:nnz-1)
        real(kind=8), intent(in)  :: vals(0:nnz-1)
        real(kind=8), intent(in)  :: r(0:ntot-1)
        real(kind=8), intent(out) :: y(0:ntot-1)
        integer :: d

        if (ntot > bgs_ntot) then          ! grow persistent buffers if needed
            if (allocated(Ay)) deallocate(Ay, rw)
            allocate(Ay(0:ntot-1), rw(0:ntot-1))
            bgs_ntot = ntot
        endif
        y = 0.d0
        ! forward sweep: u then v
        call apply_line_solve(D0, r, y)                              ! y_u = A00^-1 r_u
        call sub_matvec(rowPtr, cols, vals, nnz, ntot, 1, 0, y, Ay)  ! Ay_v = A10 y_u
        rw = r
        do d = 1, ntot-1, 2                                          ! v-dofs (odd)
            rw(d) = r(d) - Ay(d)
        end do
        call apply_line_solve(D1, rw, y)                            ! y_v = A11^-1(r_v - A10 y_u)
        if (bgs_symmetric) then
            ! backward sweep: re-solve u accounting for A01 y_v (stronger PC)
            call sub_matvec(rowPtr, cols, vals, nnz, ntot, 0, 1, y, Ay)  ! Ay_u = A01 y_v
            rw = r
            do d = 0, ntot-1, 2                                      ! u-dofs (even)
                rw(d) = r(d) - Ay(d)
            end do
            call apply_line_solve(D0, rw, y)                        ! y_u = A00^-1(r_u - A01 y_v)
        endif
    end subroutine apply_block_gs

    ! -----------------------------------------------------------------
    ! Once per level, verify the line extraction + Thomas solve on the
    ! REAL assembled matrix by a round trip:  A_block * (A_block^-1 r) == r
    ! for each field.  Prints diagnostics (path count, line lengths,
    ! smallest diagonal, round-trip error, cycle flag).  Read-only on the
    ! matrix; safe under the MPI solver server (runs on the main rank).
    ! -----------------------------------------------------------------
    subroutine validate_line_solve(rowPtr, cols, vals, nnz, ntot, level)
        integer,      intent(in) :: nnz, ntot, level
        integer,      intent(in) :: rowPtr(0:ntot), cols(0:nnz-1)
        real(kind=8), intent(in) :: vals(0:nnz-1)

        integer, parameter :: MAXLEV = 30
        logical, save :: done(MAXLEV) = .false.
        type(line_decomp_t) :: D
        real(kind=8), allocatable :: r0(:), y(:), chk(:)
        integer :: par, r, np, plen, minlen, maxlen
        real(kind=8) :: err, mindiag

        if (.not. tridiag_validate) return
        if (level >= 1 .and. level <= MAXLEV) then
            if (done(level)) return
            done(level) = .true.
        endif

        allocate(r0(0:ntot-1), y(0:ntot-1), chk(0:ntot-1))
        write(*,'(a,i0,a,i0)') ' [tridiag] validating level ', level, '  ntot=', ntot

        do par = 0, 1
            call build_line_decomp(rowPtr, cols, vals, nnz, ntot, par, D)
            if (D%nfield == 0) then
                write(*,*) '   field', par, ': empty (no bouss cells)'
                call free_line_decomp(D)
                cycle
            endif

            ! deterministic rhs on this field's dofs
            r0 = 0.d0
            do r = par, ntot-1, 2
                r0(r) = 1.d0 + 0.5d0*sin(real(r+1,8))
            end do
            y = 0.d0
            call apply_line_solve(D, r0, y)                      ! y = A_block^-1 r0
            call sub_matvec(rowPtr, cols, vals, nnz, ntot, par, par, y, chk) ! chk = A_block y

            err = 0.d0
            do r = par, ntot-1, 2
                err = max(err, abs(chk(r) - r0(r)))
            end do

            minlen = huge(1); maxlen = 0
            do np = 1, D%npaths
                plen = D%pathStart(np+1) - D%pathStart(np)
                minlen = min(minlen, plen); maxlen = max(maxlen, plen)
            end do
            mindiag = huge(1.d0)
            do r = 1, D%nfield
                mindiag = min(mindiag, abs(D%diag(r)))
            end do

            write(*,*) '   field', par, ' npaths', D%npaths, ' nfield', D%nfield, &
                       ' linelen', minlen, maxlen, ' mindiag', mindiag, &
                       ' roundtrip_err', err, ' cycle', D%has_cycle
            call free_line_decomp(D)
        end do
        deallocate(r0, y, chk)
    end subroutine validate_line_solve

end module bouss_tridiag_module
