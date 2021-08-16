program diagonalize
    use, intrinsic :: iso_fortran_env
    implicit none
    INTEGER :: m, samples, iterations
    parameter (samples = 100)
    parameter (iterations = 5)
    DOUBLE PRECISION, allocatable :: matrix(:,:), eigvecs(:,:)
    DOUBLE PRECISION, allocatable :: eigvals(:)
    INTEGER :: i,j,k,l
    REAL :: startT, endT  !! timing parameters
    DOUBLE PRECISION, allocatable :: execT(:)     !! timing vector
    DOUBLE PRECISION :: rn       !! random number

    allocate(execT(iterations))
    execT = 0.0
    m = 200

    call random_seed()
    do i = 1, iterations
        do j = 1, samples
            allocate(matrix(m,m))
            allocate(eigvecs(m,m))
            allocate(eigvals(m))
            matrix = 0.0
            do k = 1, m
                do l = k, m
                    call random_number(rn)
                    matrix(l,k) = rn
                enddo
            enddo
            call cpu_time(startT)
            call hermitean_diagonalization(matrix,m,eigvals,eigvecs)
            call cpu_time(endT)
            execT(i) = execT(i) + (endT - startT)
            deallocate(matrix)
            deallocate(eigvals)
            deallocate(eigvecs)
            write(*,*) 'm = ', m, 'iteration: ', j
        enddo
        m = m*2
    enddo
    write(*,*) execT

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
