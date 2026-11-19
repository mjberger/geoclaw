! =====================================================================
! bouss_gmres_module
!
! Self-contained, OpenMP-threaded linear solver to replace the PETSc KSP
! for the SGN Boussinesq system.  Everything in the Krylov loop is threaded:
!   - CSR matvec  (parallel over rows)
!   - vector ops  (dot/axpy/scale as OpenMP reductions/loops)
!   - preconditioner: the block Gauss-Seidel OpenMP Thomas apply
!     (apply_block_gs in bouss_tridiag_module)
!
! Right-preconditioned restarted GMRES(m).  Right preconditioning makes the
! Arnoldi residual equal the TRUE residual ||b - A x||, so the stopping test
! matches PETSc's default (||r_k|| <= rtol*||r_0||, with x_0 = 0 => ||r_0||=||b||).
!
! PETSc-free: depends only on OpenMP and bouss_tridiag_module, so it
! unit-tests with plain gfortran.
! =====================================================================
module bouss_gmres_module

    use bouss_tridiag_module, only: line_decomp_t, apply_block_gs
    implicit none
    private
    public :: csr_matvec, gmres_solve

    ! persistent work buffers for gmres_solve, grown as needed and reused across
    ! solves (avoids re-allocating the ~80 MB Krylov basis V every solve).
    ! gmres_solve is entered on the master thread only, so save is thread-safe.
    real(kind=8), allocatable, save :: V(:,:), H(:,:)
    real(kind=8), allocatable, save :: cs(:), sn(:), g(:), yy(:)
    real(kind=8), allocatable, save :: w(:), z(:), zsum(:), reshist(:)
    integer, save :: wk_ntot = -1, wk_m = -1

    ! Orthogonalization of the new Krylov vector w against V(:,1:j).  Selected
    ! at run time with the environment variable BOUSS_GMRES_ORTH (read once):
    !   1 = modified Gram-Schmidt, each dot/axpy its own OpenMP region
    !       (2j fork/joins per iteration; the original implementation)
    !   2 = modified Gram-Schmidt, serial (no OpenMP in the orthogonalization)
    !   3 = classical Gram-Schmidt in two threaded block passes (h = V^T w, then
    !       w = w - V h), re-orthogonalized once when the update removed more
    !       than a factor 1/sqrt(2) of ||w|| (the DGKS "twice is enough" test,
    !       also PETSc's default).  2 to 4 regions per iteration instead of 2j.
    integer, public, save :: gmres_orth = 3
    integer, public, save :: gmres_nreorth = 0    ! CGS re-orthogonalizations so far
    ! Mode 3 re-orthogonalization test: refine when ||w_after|| < tol*||w_before||.
    ! tol = 1/sqrt(2) is the DGKS criterion (also PETSc's): it guarantees
    ! eps-level orthogonality, but in GMRES on a well-preconditioned system the
    ! new vector lies mostly in the old span, so it fires on nearly every
    ! iteration and the cost doubles (CGS2).  GMRES only needs the loss of
    ! orthogonality, about eps/tol, to stay far below rtol, so the default is
    ! 1e-3 (margin ~1e4 at rtol 1e-9).  Environment variable BOUSS_CGS_REORTH
    ! overrides (set 0.7071 for the DGKS behaviour).
    real(kind=8), public, save :: cgs_reorth_tol = 1.d-3
    ! doubleCheck (default on): a convergence claim made inside a cycle from the
    ! Givens-rotated residual estimate is verified with the true residual
    ! ||b - A x|| (one extra matvec per solve, ~2-4%) before gmres_solve
    ! returns; if the true residual is still above rtol the iteration simply
    ! continues with another restart cycle.  This makes the stopping rule
    ! independent of how well the Krylov basis kept its orthogonality.  Set it
    ! .false. here, or via environment variable BOUSS_GMRES_DOUBLECHECK=0, to
    ! trust the estimate alone (the behaviour before Nov 2026).
    logical, public, save :: doubleCheck = .true.
    integer, public, save :: gmres_ncheckfail = 0     ! verifications that failed
    real(kind=8), public, save :: gmres_max_est_err = 0.d0, gmres_max_true_resid = 0.d0
    logical, save :: orth_init = .false.
    real(kind=8), allocatable, save :: hpart(:,:) ! (m+1, nthreads) partial sums, mode 3
    integer, save :: hp_m = -1, hp_nth = -1

contains

    ! y = A x   (compressed CSR, 0-based), threaded over rows
    subroutine csr_matvec(rowPtr, cols, vals, nnz, ntot, x, y)
        integer,      intent(in)  :: nnz, ntot, rowPtr(0:ntot), cols(0:nnz-1)
        real(kind=8), intent(in)  :: vals(0:nnz-1), x(0:ntot-1)
        real(kind=8), intent(out) :: y(0:ntot-1)
        integer :: r, p
        real(kind=8) :: s
        !$omp parallel do schedule(static) default(shared) private(r,p,s)
        do r = 0, ntot-1
            s = 0.d0
            do p = rowPtr(r), rowPtr(r+1)-1
                s = s + vals(p)*x(cols(p))
            end do
            y(r) = s
        end do
        !$omp end parallel do
    end subroutine csr_matvec

    real(kind=8) function omp_dot(n, x, y) result(d)
        integer,      intent(in) :: n
        real(kind=8), intent(in) :: x(0:n-1), y(0:n-1)
        integer :: i
        d = 0.d0
        !$omp parallel do reduction(+:d) default(shared) private(i)
        do i = 0, n-1
            d = d + x(i)*y(i)
        end do
        !$omp end parallel do
    end function omp_dot

    subroutine omp_axpy(n, a, x, y)      ! y = y + a*x
        integer,      intent(in)    :: n
        real(kind=8), intent(in)    :: a, x(0:n-1)
        real(kind=8), intent(inout) :: y(0:n-1)
        integer :: i
        !$omp parallel do default(shared) private(i)
        do i = 0, n-1
            y(i) = y(i) + a*x(i)
        end do
        !$omp end parallel do
    end subroutine omp_axpy

    subroutine omp_scalecopy(n, a, x, y)  ! y = a*x
        integer,      intent(in)  :: n
        real(kind=8), intent(in)  :: a, x(0:n-1)
        real(kind=8), intent(out) :: y(0:n-1)
        integer :: i
        !$omp parallel do default(shared) private(i)
        do i = 0, n-1
            y(i) = a*x(i)
        end do
        !$omp end parallel do
    end subroutine omp_scalecopy

    ! apply the preconditioner: z = M^{-1} r  (block G-S, or identity)
    subroutine apply_pc(precondition, D0, D1, rowPtr, cols, vals, nnz, ntot, r, z)
        logical,             intent(in)  :: precondition
        type(line_decomp_t), intent(in)  :: D0, D1
        integer,             intent(in)  :: nnz, ntot, rowPtr(0:ntot), cols(0:nnz-1)
        real(kind=8),        intent(in)  :: vals(0:nnz-1), r(0:ntot-1)
        real(kind=8),        intent(out) :: z(0:ntot-1)
        if (precondition) then
            call apply_block_gs(D0, D1, rowPtr, cols, vals, nnz, ntot, r, z)
        else
            z = r
        endif
    end subroutine apply_pc

    ! -----------------------------------------------------------------
    ! Right-preconditioned restarted GMRES(m), zero initial guess.
    !   A given as compressed CSR (rowPtr,cols,vals); the preconditioner is
    !   block Gauss-Seidel using the line decompositions D0 (u), D1 (v)
    !   (set precondition=.false. for unpreconditioned GMRES).
    !   Converges when ||b - A x||_2 <= rtol * ||b||_2.
    !   iters = total inner iterations, resid = final relative residual,
    !   info  = 0 converged, 1 hit maxcycles*m without converging.
    ! -----------------------------------------------------------------
    subroutine gmres_solve(rowPtr, cols, vals, nnz, ntot, D0, D1, precondition, &
                           b, x, rtol, m, maxcycles, iters, resid, info)
        integer,      intent(in)  :: nnz, ntot, rowPtr(0:ntot), cols(0:nnz-1)
        real(kind=8), intent(in)  :: vals(0:nnz-1)
        type(line_decomp_t), intent(in) :: D0, D1
        logical,      intent(in)  :: precondition
        real(kind=8), intent(in)  :: b(0:ntot-1)
        real(kind=8), intent(out) :: x(0:ntot-1)
        real(kind=8), intent(in)  :: rtol
        integer,      intent(in)  :: m, maxcycles
        integer,      intent(out) :: iters, info
        real(kind=8), intent(out) :: resid

        real(kind=8) :: beta, bnorm, denom, tmp, s
        logical, save :: hist_printed = .false.
        integer :: i, j, k, cyc, jj, r, nth, istat, ios
        character(len=16) :: envstr
        logical :: est_conv
        integer :: omp_get_max_threads      ! external, called only on !$ lines

        if (.not. orth_init) then           ! pick the orthogonalization once
            call get_environment_variable('BOUSS_GMRES_ORTH', envstr, status=istat)
            if (istat == 0 .and. len_trim(envstr) > 0) then
                read(envstr, *, iostat=ios) gmres_orth
                if (ios /= 0) gmres_orth = 3
            endif
            if (gmres_orth < 1 .or. gmres_orth > 3) gmres_orth = 3
            call get_environment_variable('BOUSS_CGS_REORTH', envstr, status=istat)
            if (istat == 0 .and. len_trim(envstr) > 0) then
                read(envstr, *, iostat=ios) cgs_reorth_tol
                if (ios /= 0) cgs_reorth_tol = 1.d-3
            endif
            call get_environment_variable('BOUSS_GMRES_DOUBLECHECK', envstr, status=istat)
            if (istat == 0 .and. len_trim(envstr) > 0) doubleCheck = (envstr(1:1) /= '0')
            write(*,'(a,i0,a)') ' GMRES orthogonalization mode ', gmres_orth, &
                  '  (1 = MGS threaded, 2 = MGS serial, 3 = CGS + re-orth test)'
            if (gmres_orth == 3) write(*,'(a,es9.2)') &
                  ' GMRES CGS re-orthogonalization threshold ', cgs_reorth_tol
            write(*,'(a,l1)') ' GMRES doubleCheck (true residual verified at convergence): ', &
                  doubleCheck
            orth_init = .true.
        endif

        if (ntot > wk_ntot .or. m > wk_m) then   ! grow persistent buffers
            if (allocated(V)) deallocate(V,H,cs,sn,g,yy,w,z,zsum,reshist)
            allocate(V(0:ntot-1, m+1), H(m+1, m))
            allocate(cs(m), sn(m), g(m+1), yy(m))
            allocate(w(0:ntot-1), z(0:ntot-1), zsum(0:ntot-1))
            allocate(reshist(0:m))
            wk_ntot = ntot; wk_m = m
        endif
        nth = 1
        !$ nth = omp_get_max_threads()
        if (m+1 > hp_m .or. nth > hp_nth) then  ! per-thread partials for mode 3
            if (allocated(hpart)) deallocate(hpart)
            allocate(hpart(m+1, nth))
            hp_m = m+1; hp_nth = nth
        endif
        reshist = 0.d0

        x = 0.d0
        iters = 0
        info  = 1
        bnorm = sqrt(omp_dot(ntot, b, b))
        if (bnorm == 0.d0) then
            resid = 0.d0; info = 0
            return
        endif

        cyc = 0
        est_conv = .false.
        do
            ! true residual r = b - A x  (x0=0 the first time).  With doubleCheck
            ! the loop also comes back here after the in-cycle estimate reached
            ! rtol, so the claim is verified before returning.
            call csr_matvec(rowPtr, cols, vals, nnz, ntot, x, w)
            call omp_axpy(ntot, -1.d0, b, w)      ! w = A x - b
            beta = sqrt(omp_dot(ntot, w, w))      ! ||A x - b|| = ||r||
            if (est_conv) then                    ! verifying the estimate in `resid`
                gmres_max_est_err    = max(gmres_max_est_err, abs(beta/bnorm - resid))
                gmres_max_true_resid = max(gmres_max_true_resid, beta/bnorm)
                if (beta > rtol*bnorm) gmres_ncheckfail = gmres_ncheckfail + 1
                est_conv = .false.
            endif
            resid = beta / bnorm
            reshist(0) = resid
            if (beta <= rtol*bnorm) then; info = 0; exit; endif
            if (cyc >= maxcycles) exit            ! out of cycles; info stays 1
            cyc = cyc + 1
            call omp_scalecopy(ntot, -1.d0/beta, w, V(:,1))   ! v1 = (b - A x)/beta
            g = 0.d0; g(1) = beta

            k = 0
            do j = 1, m
                k = j
                call apply_pc(precondition, D0, D1, rowPtr, cols, vals, nnz, &
                              ntot, V(:,j), z)                        ! z = M^-1 v_j
                call csr_matvec(rowPtr, cols, vals, nnz, ntot, z, w)  ! w = A z
                ! orthogonalize w against V(:,1:j); sets H(1:j,j), H(j+1,j)=||w||
                select case (gmres_orth)
                case (2)        ! modified Gram-Schmidt, serial
                    do i = 1, j
                        s = 0.d0
                        do r = 0, ntot-1
                            s = s + w(r)*V(r,i)
                        end do
                        H(i,j) = s
                        do r = 0, ntot-1
                            w(r) = w(r) - s*V(r,i)
                        end do
                    end do
                    s = 0.d0
                    do r = 0, ntot-1
                        s = s + w(r)*w(r)
                    end do
                    H(j+1,j) = sqrt(s)
                case (3)        ! classical Gram-Schmidt, block passes + re-orth test
                    call cgs_orth(ntot, j)
                case default    ! modified Gram-Schmidt, one OpenMP region per dot/axpy
                    do i = 1, j
                        H(i,j) = omp_dot(ntot, w, V(:,i))
                        call omp_axpy(ntot, -H(i,j), V(:,i), w)
                    end do
                    H(j+1,j) = sqrt(omp_dot(ntot, w, w))
                end select
                if (H(j+1,j) > 0.d0) then
                    call omp_scalecopy(ntot, 1.d0/H(j+1,j), w, V(:,j+1))
                else
                    V(:,j+1) = 0.d0
                endif
                ! apply previous Givens rotations to new column
                do i = 1, j-1
                    tmp      =  cs(i)*H(i,j) + sn(i)*H(i+1,j)
                    H(i+1,j) = -sn(i)*H(i,j) + cs(i)*H(i+1,j)
                    H(i,j)   =  tmp
                end do
                ! new Givens rotation to eliminate H(j+1,j)
                denom = sqrt(H(j,j)**2 + H(j+1,j)**2)
                cs(j) = H(j,j)/denom
                sn(j) = H(j+1,j)/denom
                H(j,j)   = cs(j)*H(j,j) + sn(j)*H(j+1,j)
                H(j+1,j) = 0.d0
                tmp    =  cs(j)*g(j)
                g(j+1) = -sn(j)*g(j)
                g(j)   =  tmp
                resid  = abs(g(j+1)) / bnorm
                reshist(j) = resid
                iters  = iters + 1
                if (abs(g(j+1)) <= rtol*bnorm) exit
            end do

            ! solve upper-triangular H(1:k,1:k) yy = g(1:k)
            do i = k, 1, -1
                tmp = g(i)
                do j = i+1, k
                    tmp = tmp - H(i,j)*yy(j)
                end do
                yy(i) = tmp / H(i,i)
            end do
            ! x += M^{-1} (V(:,1:k) yy)   [right preconditioning]
            zsum = 0.d0
            do i = 1, k
                call omp_axpy(ntot, yy(i), V(:,i), zsum)
            end do
            call apply_pc(precondition, D0, D1, rowPtr, cols, vals, nnz, ntot, zsum, z)
            call omp_axpy(ntot, 1.d0, z, x)

            if (resid <= rtol) then               ! Givens estimate says converged
                if (.not. doubleCheck) then
                    info = 0; exit
                endif
                est_conv = .true.                 ! go verify with the true residual
            else
                ! only prints when a full cycle finished without converging
                write(*,'(a,i3,a,i5,a,es11.3)') " GMRES restart cycle ",cyc,   &
                      "  total its ",iters,"  rel resid ",resid
                if (.not. hist_printed) then
                    hist_printed = .true.
                    write(*,'(a,es12.4,a,i0)') "  [gmres stall] ||b|| = ",bnorm, &
                          "   ntot = ",ntot
                    write(*,'(a)') "  [gmres stall] intra-cycle relative residual:"
                    do jj = 0, k
                        if (jj <= 6 .or. mod(jj,10) == 0)                        &
                          write(*,'(a,i4,a,es12.4)') "      iter ",jj,":  ",reshist(jj)
                    end do
                endif
            endif
        end do
    end subroutine gmres_solve

    ! -----------------------------------------------------------------
    ! [mode 3] Classical Gram-Schmidt of the module vector w against
    ! V(:,1:j) in two block passes: all j dot products in one threaded sweep,
    ! then all j updates in one threaded sweep (each thread keeps its chunk of
    ! w in cache and streams each column of V once).  DGKS test: if the update
    ! removed more than a factor 1/sqrt(2) of ||w||, the coefficients may have
    ! lost accuracy to cancellation, so do it a second time ("twice is
    ! enough").  On exit H(1:j,j) holds the coefficients and H(j+1,j) = ||w||.
    ! -----------------------------------------------------------------
    subroutine cgs_orth(ntot, j)
        integer, intent(in) :: ntot, j
        real(kind=8) :: hloc(j), wn0, wn1

        call block_dots(ntot, j, hloc, wn0)       ! hloc = V^T w,  wn0 = ||w||^2 before
        call block_update(ntot, j, hloc, wn1)     ! w = w - V hloc, wn1 = ||w||^2 after
        H(1:j,j) = hloc
        if (wn1 < cgs_reorth_tol**2 * wn0) then
            gmres_nreorth = gmres_nreorth + 1
            call block_dots(ntot, j, hloc, wn0)
            call block_update(ntot, j, hloc, wn1)
            H(1:j,j) = H(1:j,j) + hloc
        endif
        H(j+1,j) = sqrt(wn1)
    end subroutine cgs_orth

    ! rows 0:ntot-1 split into nth contiguous chunks; chunk tid is r0:r1
    pure subroutine chunk_bounds(ntot, nth, tid, r0, r1)
        integer, intent(in)  :: ntot, nth, tid
        integer, intent(out) :: r0, r1
        r0 = int((int(ntot,8)*tid)/nth)
        r1 = int((int(ntot,8)*(tid+1))/nth) - 1
    end subroutine chunk_bounds

    ! h(i) = V(:,i).w for i = 1..j and wn = w.w, in ONE threaded sweep.
    ! Per-thread partial sums go to hpart(:,tid+1) and are combined in thread
    ! order afterwards, so the result is reproducible for a fixed thread count.
    subroutine block_dots(ntot, j, h, wn)
        integer,      intent(in)  :: ntot, j
        real(kind=8), intent(out) :: h(j), wn
        integer :: i, r, r0, r1, nth, tid, t, nteam
        real(kind=8) :: s
        integer :: omp_get_num_threads, omp_get_thread_num   ! external, !$ lines only

        nteam = 1
        !$omp parallel default(shared) private(i,r,r0,r1,nth,tid,s)
        nth = 1; tid = 0
        !$ nth = omp_get_num_threads()
        !$ tid = omp_get_thread_num()
        if (tid == 0) nteam = nth
        call chunk_bounds(ntot, nth, tid, r0, r1)
        do i = 1, j
            s = 0.d0
            do r = r0, r1
                s = s + V(r,i)*w(r)
            end do
            hpart(i, tid+1) = s
        end do
        s = 0.d0
        do r = r0, r1
            s = s + w(r)*w(r)
        end do
        hpart(j+1, tid+1) = s
        !$omp end parallel

        h = 0.d0; wn = 0.d0
        do t = 1, nteam
            h(1:j) = h(1:j) + hpart(1:j, t)
            wn = wn + hpart(j+1, t)
        end do
    end subroutine block_dots

    ! w = w - sum_i h(i) V(:,i)  and  wn = ||w||^2 afterwards, ONE threaded sweep
    subroutine block_update(ntot, j, h, wn)
        integer,      intent(in)  :: ntot, j
        real(kind=8), intent(in)  :: h(j)
        real(kind=8), intent(out) :: wn
        integer :: i, r, r0, r1, nth, tid, t, nteam
        real(kind=8) :: s, hi
        integer :: omp_get_num_threads, omp_get_thread_num   ! external, !$ lines only

        nteam = 1
        !$omp parallel default(shared) private(i,r,r0,r1,nth,tid,s,hi)
        nth = 1; tid = 0
        !$ nth = omp_get_num_threads()
        !$ tid = omp_get_thread_num()
        if (tid == 0) nteam = nth
        call chunk_bounds(ntot, nth, tid, r0, r1)
        do i = 1, j
            hi = h(i)
            do r = r0, r1
                w(r) = w(r) - hi*V(r,i)
            end do
        end do
        s = 0.d0
        do r = r0, r1
            s = s + w(r)*w(r)
        end do
        hpart(j+1, tid+1) = s
        !$omp end parallel

        wn = 0.d0
        do t = 1, nteam
            wn = wn + hpart(j+1, t)
        end do
    end subroutine block_update

end module bouss_gmres_module
