! Biconjugate gradient stabilized method
! COPMPILATION: f90 linalg.f90 bicg.f90 -o bicg.exe -xopenmp
! RUN         : ./bicg.exe Ndim NumProc -key
!                               Ndim - matrix dimension; NumProc - number of openmp threads
!                               key: -3diag - test run for 3 diag matrix
!                                        -identity - test run for E matrix

program bicg
use linalg
use omp_lib
implicit none

character(len=32) :: arg, mode
integer :: N = 10, NUMPROC = 1
integer :: i,j,k,it 
real*8, dimension(:,:), allocatable :: A
real*8, dimension(:), allocatable :: b, x, s, t, r, p, v
real*8 :: alpha, beta, rho0, rho1, omega, bb, rr

real*8, parameter :: eps = 1.0E-5
real*8 :: t1, t2, t3

! command arguments
do
    call get_command_argument(i,arg)
    if (len_trim(arg) == 0) exit
    if (i == 1) read(arg(:),*) N ! dimension of matrix
    if (i == 2) read(arg(:),*) NUMPROC  ! num of processes
    if (i == 3) mode = trim(arg)
    i = i + 1
enddo

print '("Dimension: ",i7)',N
print '("Num of processes: ",i3)', NUMPROC
print '("Accuracy: ", ES10.2)', eps

t1 = omp_get_wtime()

allocate(A(n,n))
A = 0.
allocate(b(n))
b = 0.
 
 call omp_set_num_threads(NUMPROC)
! check args
if (trim(mode) == '-identity') mode = '-test'
select case (trim(mode))
    case('-test')
       print *, "Test mode is enable: A=E, b=x."
       !$OMP PARALLEL DO PRIVATE(i), SHARED(A,b)
       do i = 1, n
            A(i,i) = 1.
            b(i) = i * 1.
       enddo
    case('-3diag')
       print *, "Test mode is enable: A is 3-diag matrix"
       !$OMP PARALLEL DO PRIVATE(i), SHARED(A,b)
       do i = 2, n-1
            A(i-1:i+1,i) = (/1.0,-3.0,1.0/)
            b(i) = -5.
         enddo
       A(1:2,1) = (/-3.0,1.0/); A(N-1:N,N) = (/1.0,-3.0/)
       b(1) = -10.; b(n) = -10.
        case default
          print *, "Try -identity or -3diag flags. Run for bad nonsymmetric matrix."
          A = 1.333
          A(1,:n/3) = 100.
          A(n,:n/3) = -10.
          b = n * 1.333
          b(1) = 0.1
          b(n) = -133.901
end select

t2 = omp_get_wtime()

allocate(s(n))
allocate(t(n))
allocate(x(n))
allocate(r(n))
allocate(p(n))
allocate(v(n))
x = 0.
p = 0.
v = 0.
r = b
rho0 = 1.
alpha = 1.
omega = 1.
bb = norm2(b)
rr = bb

! BiCG: let r0 = b
it = 0
do while (rr / bb > eps)
        rho1 = dotprod(b,r)
        beta = rho1 * alpha / (rho0 * omega)
        p = r + beta * (p - omega * v)
        v = multmv(A,p)
        alpha = rho1 / dotprod(b, v)
        s = r - alpha * v
        t = multmv(A,s)
        omega = dotprod(t,s) / dotprod(t,t)
        x = x + alpha * p + omega * s
        r = s - omega * t
        rr = norm2(r)
        rho0 = rho1
    it = it + 1
    if (mod(it,250) == 0) print"('DEBUG ||r||2/||b||2 = ',ES10.2, ' on ',i5 ' iters.')", rr / bb, it
enddo

t3 = omp_get_wtime()
print "('||r||2/||b||2 = ',ES10.2, ' on ',i5 ' iters.')", rr / bb, it
print '("Method time ",f10.5,"s for ",i7,"-size on ", i3 " threads.")', t3 - t2, N, NUMPROC
print '("Initial time, s: ",f10.5)',t2 - t1
print '("CHECK of solve (method): ",L)',check_solve(a,x,b,eps)
if(N <= 10)print *,"x: ", x

deallocate(A)
deallocate(b)
deallocate(x)
deallocate(r)
deallocate(v)
deallocate(p)
deallocate(s)
deallocate(t)
end program bicg
