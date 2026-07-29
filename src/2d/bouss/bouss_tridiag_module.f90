! =====================================================================
! bouss_tridiag_module
!
! Item 1-2 of the parallel PCSHELL preconditioner for the SGN Boussinesq
! linear system.  The 2x2-block matrix (interleaved u0,v0,u1,v1,...) has
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
    public :: validate_line_solve

    ! set .false. to disable the (once-per-level) round-trip validation check
    logical, public :: tridiag_validate = .true.

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
        logical :: has_cycle = .false.        ! true if a periodic (cyclic) line was found
    end type line_decomp_t

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
        integer :: np, slot, cur, prev, nxt, startf, ncount

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

        ! ---- build symmetric-closure adjacency among field dofs ----
        allocate(nbr(2,nfield)); nbr = -1
        allocate(deg(nfield));   deg = 0

        do f = 1, nfield
            r = f2dof(f)
            do p = rowPtr(r), rowPtr(r+1)-1
                c = cols(p)
                if (c == r) cycle                 ! diagonal, not an edge
                if (c < 0 .or. c > ntot-1) cycle
                if (mod(c,2) /= parity) cycle      ! other field (cross-coupling), skip
                if (vals(p) == 0.d0) cycle         ! explicit zero, not an edge
                cf = dof2f(c)
                call add_nbr(nbr, deg, nfield, f, cf)  ! symmetric: adds both directions
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

        ! ---- extract tridiagonal coefficients along each ordered path ----
        allocate(D%sub(nfield), D%diag(nfield), D%sup(nfield))
        D%sub = 0.d0; D%diag = 0.d0; D%sup = 0.d0
        do np = 1, D%npaths
            ncount = D%pathStart(np+1) - D%pathStart(np)
            do slot = D%pathStart(np), D%pathStart(np+1)-1
                r = D%dof(slot)
                D%diag(slot) = get_entry(rowPtr, cols, vals, nnz, ntot, r, r)
                if (slot > D%pathStart(np)) then
                    prev = D%dof(slot-1)
                    D%sub(slot) = get_entry(rowPtr, cols, vals, nnz, ntot, r, prev)
                endif
                if (slot < D%pathStart(np+1)-1) then
                    nxt = D%dof(slot+1)
                    D%sup(slot) = get_entry(rowPtr, cols, vals, nnz, ntot, r, nxt)
                endif
            end do
        end do

        deallocate(dof2f, f2dof, nbr, deg, visited)
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
    subroutine apply_line_solve(D, rhs, x)
        type(line_decomp_t), intent(in)  :: D
        real(kind=8),        intent(in)  :: rhs(0:D%ntot-1)
        real(kind=8),        intent(out) :: x(0:D%ntot-1)
        integer :: p

        !$omp parallel do schedule(dynamic,16) default(shared) private(p)
        do p = 1, D%npaths
            call solve_one_path(D, p, rhs, x)
        end do
        !$omp end parallel do
    end subroutine apply_line_solve

    subroutine solve_one_path(D, p, rhs, x)
        type(line_decomp_t), intent(in)  :: D
        integer,             intent(in)  :: p
        real(kind=8),        intent(in)  :: rhs(0:D%ntot-1)
        real(kind=8),        intent(out) :: x(0:D%ntot-1)
        integer :: off, n, m
        real(kind=8) :: a(D%pathStart(p+1)-D%pathStart(p))
        real(kind=8) :: b(D%pathStart(p+1)-D%pathStart(p))
        real(kind=8) :: c(D%pathStart(p+1)-D%pathStart(p))
        real(kind=8) :: r(D%pathStart(p+1)-D%pathStart(p))
        real(kind=8) :: sol(D%pathStart(p+1)-D%pathStart(p))

        off = D%pathStart(p) - 1
        n   = D%pathStart(p+1) - D%pathStart(p)
        do m = 1, n
            a(m) = D%sub(off+m)
            b(m) = D%diag(off+m)
            c(m) = D%sup(off+m)
            r(m) = rhs(D%dof(off+m))
        end do
        call thomas(n, a, b, c, r, sol)
        do m = 1, n
            x(D%dof(off+m)) = sol(m)
        end do
    end subroutine solve_one_path

    ! Standard Thomas algorithm for a tridiagonal system.
    !   a = subdiagonal (a(1) unused), b = diagonal,
    !   c = superdiagonal (c(n) unused), d = rhs, x = solution.
    subroutine thomas(n, a, b, c, d, x)
        integer,      intent(in)  :: n
        real(kind=8), intent(in)  :: a(n), b(n), c(n), d(n)
        real(kind=8), intent(out) :: x(n)
        real(kind=8) :: cp(n), dp(n), denom
        integer :: i

        cp(1) = c(1) / b(1)
        dp(1) = d(1) / b(1)
        do i = 2, n
            denom = b(i) - a(i)*cp(i-1)
            cp(i) = c(i) / denom
            dp(i) = (d(i) - a(i)*dp(i-1)) / denom
        end do
        x(n) = dp(n)
        do i = n-1, 1, -1
            x(i) = dp(i) - cp(i)*x(i+1)
        end do
    end subroutine thomas

    subroutine free_line_decomp(D)
        type(line_decomp_t), intent(inout) :: D
        if (allocated(D%pathStart)) deallocate(D%pathStart)
        if (allocated(D%dof))       deallocate(D%dof)
        if (allocated(D%sub))       deallocate(D%sub)
        if (allocated(D%diag))      deallocate(D%diag)
        if (allocated(D%sup))       deallocate(D%sup)
        D%ntot = 0; D%parity = -1; D%nfield = 0; D%npaths = 0
        D%has_cycle = .false.
    end subroutine free_line_decomp

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
