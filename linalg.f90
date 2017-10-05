module linalg
use omp_lib
implicit none
public

contains

! Vectors multiplication
real(kind = 8) function dotprod(u, v)
real(kind = 8), dimension(:), intent(in) :: u, v
integer :: i
dotprod = 0.
if(size(u) /= size(v))return
!$OMP PARALLEL DO REDUCTION(+:dotprod)
do i=1,size(u)
    dotprod = dotprod + u(i) * v(i)
enddo
end function dotprod

! inf norm = maximum norm
real(kind = 8) function norminf(u)
real(kind = 8), dimension(:), intent(in) :: u
integer :: i
norminf = abs(u(1))
!$OMP PARALLEL DO REDUCTION(max:norminf)
do i = 2, size(u)
	if (abs(u(i)) > norminf) norminf = abs(u(i))
enddo
end function norminf

! L1 norm = Manhattan norm
real(kind = 8) function norm1(u)
real(kind = 8), dimension(:), intent(in) :: u
integer :: i
norm1 = 0.
!$OMP PARALLEL DO REDUCTION(+:norm1)
do i = 1, size(u)
    norm1 = norm1 + abs(u(i))
enddo
end function norm1

! L2 norm = Euclidian distance
real(kind = 8) function norm2(u)
real(kind = 8), dimension(:), intent(in) :: u
norm2 = sqrt(dotprod(u, u))
end function norm2

! Matrix multiplication
function multmat(A, B)
real(kind = 8), dimension(:,:), intent(in) :: A, B
real(kind = 8), dimension(size(A, dim = 1), size(B, dim = 2))  :: multmat
real(kind = 8), dimension(size(A, dim = 2)) :: row
integer :: i,j,k,m,n,p

multmat = 0.
n = size(A, dim = 1)
m = size(B, dim = 2)
p = size(A, dim = 2)
if(size(A, dim = 2) /= size(B, dim = 1))return

!$OMP PARALLEL DO PRIVATE(i,j,k,row) SHARED(A,B,multmat,m,n,p)
do i = 1, n
    row(:) = A(i,:)
    do j = 1, m
    do k = 1, p
        multmat(i,j) = multmat(i,j) + row(k) * B(k,j)
    enddo
    enddo
enddo 
end function multmat

! Multiplicate matrix on vector
function multmv(A, v)
real(kind = 8), dimension(:,:), intent(in) :: A
real(kind = 8), dimension(:), intent(in) :: v
real(kind = 8), dimension(size(A, dim = 1))  :: multmv
integer :: i,j,m,n

multmv = 0.
n = size(A, dim = 1)
m = size(v)
if(size(A, dim = 2) /= m)return

!$OMP PARALLEL DO PRIVATE(i,j) SHARED(A,v,multmv,m,n)
do i = 1, n
    do j = 1, m
        multmv(i) = multmv(i) + A(i,j) * v(j)
    enddo
enddo 
end function multmv

! Back Gauss
pure function ReverseGauss(A, b) result(x)
real(kind = 8), dimension(:,:), intent(in) :: A
real(kind = 8), dimension(:), intent(in) :: b
real(kind = 8), dimension(size(b)) :: x
integer :: i,j,n

n = size(A, dim = 1)
if(size(b) /= n)return
x = b
do i = n, -1, -1
  do j = i+1, n
	x(i) = x(i) - x(j) * A(i,j)
  enddo
  x(i) = 1. / A(i,i) * x(i)
enddo
end function ReverseGauss

! Apply Givens rotations
subroutine ApplyGivens(c,s,v)
real(kind = 8), dimension(:), intent(in) :: c, s
real(kind = 8), dimension(:), intent(inout) :: v
real(kind = 8) :: vi
integer :: i, n

n = size(c)
if((size(s) /= n) .or. (size(v) - n /= 1)) return

do i = 1, n
  vi = v(i)
  v(i)   =  c(i) * vi + s(i) * v(i+1)
  v(i+1) = -s(i) * vi + c(i) * v(i+1)
enddo
end subroutine ApplyGivens

! Check the solution of system of linear equations
logical function check_solve(A, x, b, eps)
real(kind = 8), dimension(:,:), intent(in) :: A
real(kind = 8), dimension(:), intent(in) :: x, b
real(kind = 8), intent(in), optional :: eps
real(kind = 8) :: tol = 1.e-5
real(kind = 8) :: summ
integer :: i,j,k,n
logical :: flag

if(present(eps)) tol = eps
check_solve = .true.
n = size(A, dim = 1)
if (check_solve) then
    !$OMP PARALLEL DO PRIVATE(i,j,summ) SHARED(a,x,b,check_solve,tol,n)
    do i = 1, n
        summ = b(i)
        do j = 1, n
            summ = summ - A(i,j) * x(j)
        enddo
        if (abs(summ) > tol) then
            !$OMP CRITICAL
             check_solve = .false.
            !$OMP END CRITICAL
        endif
    enddo
    !$OMP END PARALLEL DO
endif
end function check_solve

end module linalg
