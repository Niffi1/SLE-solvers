! 
! COPMPILATION: f90 linalg.f90 gmres.f90 -o gmres.exe -xopenmp
! RUN         : ./gmres.exe Ndim NumProc -key
!				Ndim - matrix dimension; NumProc - number of openmp threads
!				key: -3diag - test run for 3 diag matrix
!					 -identity - test run for E matrix

program gmres
use linalg
use omp_lib
implicit none

character(len=32) :: arg, mode
integer :: n = 8, NUMPROC = 1, m = 5
integer :: i, j, it
real(kind = 8), dimension(:,:), allocatable :: A, V, H
real(kind = 8), dimension(:), allocatable :: b, x, r, c, s, be1
real(kind = 8) :: bb, rr, denom, be1j
real(kind = 8) :: eps = 1.0e-9
real(kind = 8) :: t1

! command arguments
do
    call get_command_argument(i,arg)
    if (len_trim(arg) == 0) exit
    if (i == 1) read(arg(:), *) n ! dimension of matrix
    if (i == 2) read(arg(:), *) NUMPROC  ! num of processes
    if (i == 3) read(arg(:), *) m ! restart every max_iter iteration 
    if (i == 4) mode = trim(arg)
    i = i + 1
enddo

print '("Dimension: ",i7)', n
print '("Num of processes: ",i3)', NUMPROC
print '("Accuracy: ", ES10.2)', eps
print '("Restart every ",i7," iterations")', m

allocate(A(n,n))
A = 0.
allocate(b(n))
b = 0.
allocate(x(n))
x = 0.
 
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
	  print *, "Try -idenitity or -3diag flags. Run for bad nonsymmetric matrix."
	  do i = 1, n
		x(i) = (i-1)*9/11.
		do j = 1, n
		  A(i,j) = 1.333 * (i-1) + j-1
		enddo
		A(i,i) = (i-1) * (i-1)
	  enddo
	  print*,'EXACT SOLUTION: ',x
	  do i = 1, n
		b(i) = 0.
		do j = 1, n
		  b(i) = b(i) + A(i,j) * x(j)
		enddo
	  enddo
end select

t1 = omp_get_wtime()
allocate(r(n))
allocate(s(m))
allocate(c(m))
allocate(be1(m+1))
allocate(V(n,m+1))
allocate(H(m+1,m))

! GMRES:
! init
x = 0.
r = b
bb = norm2(b)
if (bb == 0.) bb = 1.
rr = bb
eps = bb * eps
it = 0

! main cycle
do while (rr > eps)
	V = 0.
	H = 0.
	c = 0.
	s = 0.
	be1 = 0.
	be1(1) = rr
	
	! Arnoldi iteration
    V(:, 1) = 1. / rr * r
    do j = 1, m
  	  V(:,j+1) = multmv(A, V(:,j))
  	  do i = 1, j
		H(i,j) = dotprod(V(:,i),V(:,j+1))
		V(:,j+1) = V(:,j+1) - H(i,j) * V(:,i)
  	  enddo
	  H(j+1,j) = norm2(V(:,j+1))
  	  if(abs(H(j+1,j)) < 1.e-9) exit ! happy solution x=x(m-1)
  	  V(:,j+1) = V(:,j+1) * (1. / H(j+1,j))

  	  ! Givens rotations
  	  if (j > 1) call ApplyGivens(c(:j-1),s(:j-1),H(:j,j))
  	  
  	  denom = 1. / sqrt(H(j,j) * H(j,j) + H(j+1,j) * H(j+1,j))
  	  c(j) = H(j,j) * denom
  	  s(j) = H(j+1,j) * denom 
  	  be1j = be1(j)
  	  be1(j) = c(j) * be1j + s(j) * be1(j+1)
  	  be1(j+1) = -s(j) * be1j + c(j) * be1(j+1)
  	  H(j,j) = c(j) * H(j,j) + s(j) * H(j+1,j)
  	  H(j+1,j) = 0.
	enddo
  
    rr = abs(be1(m+1))
    c = ReverseGauss(H(:m,:),be1(:m))
	x = x + multmv(V(:,:m), c)
    r = b - multmv(A,x)
    it = it + m
    if (mod(it,100) == 0) print"('DEBUG ||r||2/||b||2 = ',ES10.2, ' on ',i5 ' iters.')", rr / bb, it
enddo

print "('||r||2/||b||2 = ',ES10.2, ' on ',i5 ' iters.')", rr / bb, it
print '("Method time ",f10.5,"s for ",i7,"-size on ", i3 " threads.")', omp_get_wtime() - t1, N, NUMPROC
print '("CHECK of solve (method): ",L)',check_solve(A,x,b,eps)
if(n <= 10)print *,"x: ", x

deallocate(A)
deallocate(V)
deallocate(H)
deallocate(be1)
deallocate(c)
deallocate(s)
deallocate(b)
deallocate(x)
deallocate(r)
end program gmres
