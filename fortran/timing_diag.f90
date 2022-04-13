program diagonalize
    use, intrinsic :: iso_fortran_env
    implicit none

    ! control/input variables
    integer :: m                ! order of matrix
    integer :: iterations       ! number of iterations, m=m*2 at each
    integer :: samples          ! number of samples at each iteration
    parameter (samples = 100)
    parameter (iterations = 5)

    ! linear algebra variables
    double precision, allocatable :: matrix(:,:), eigvecs(:,:) 
    double precision, allocatable :: eigvals(:)

    ! helpers
    integer :: i,j,k,l
    double precision :: startT, endT            ! trigger stopwatch points
    double precision, allocatable :: t(:,:)     ! timing matrix

    allocate(t(iterations,samples))
    t = 0.d0
    m = 200

    call random_seed()
    do i = 1, iterations
        allocate(matrix(m,m))
        allocate(eigvecs(m,m))
        allocate(eigvals(m))

        do j = 1, samples
            call random_number(matrix)

            call cpu_time(startT)
            call hermitean_diagonalization(matrix,m,eigvals,eigvecs)
            call cpu_time(endT)

            t(i,j) = endT - startT
        enddo

        m = m*2
        deallocate(matrix, eigvals, eigvecs)

    enddo

    open(unit=12, file="diag-timing.dat", action="write", status="new")
    do i = 1, samples
        write(12,*) t(:,i)
    enddo

    deallocate(t)
end program diagonalize

subroutine hermitean_diagonalization(matrix,ord,eigvals,eigvecs)
    IMPLICIT NONE
    INTEGER, intent(in) :: ord
    DOUBLE PRECISION, intent(in) :: matrix(ord,ord)
    DOUBLE PRECISION, intent(out) :: eigvals(ord)
    DOUBLE PRECISION, intent(out) :: eigvecs(ord,ord)
    INTEGER LWMAX
    PARAMETER (LWMAX=10000)
    INTEGER INFO, LWORK
    DOUBLE PRECISION, allocatable :: WORK(:)

    allocate(WORK(LWMAX))
    eigvecs = matrix
    LWORK=-1
    call DSYEV('Vectors', 'Upper', ord, eigvecs, ord, eigvals, WORK, LWORK, INFO)
    LWORK = min(LWMAX, int(WORK(1)))
    !write(*,*) LWORK
    deallocate(WORK)
    allocate(WORK(LWORK))
    call DSYEV('V', 'U', ord, eigvecs, ord, eigvals, WORK, LWORK, INFO)
    if(INFO == 0) then
        !write(*,*) 'Diagonalization performed.'
    else
        write(*,*) 'Diagonalization failed.'
        return
    endif

end subroutine hermitean_diagonalization
