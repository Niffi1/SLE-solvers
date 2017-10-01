module linalg
use omp_lib
implicit none
public

contains

real*8 function dotprod(u, v)
real*8, dimension(:), intent(in) :: u, v
integer :: i
dotprod = 0.
if(size(u) /= size(v))return
!$OMP PARALLEL DO REDUCTION(+:dotprod)
do i=1,size(u)
    dotprod = dotprod + u(i) * v(i)
enddo
end function dotprod

real*8 function norm2(u)
real*8, dimension(:), intent(in) :: u
norm2 = sqrt(dotprod(u, u))
end function norm2

function multmat(A, B)
real*8, dimension(:,:), intent(in) :: A, B
real*8, dimension(size(A, dim = 1), size(B, dim = 2))  :: multmat
real*8, dimension(size(A, dim = 2)) :: row
integer :: i,j,k,m,n,p
'
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

function multmv(A, v)
real*8, dimension(:,:), intent(in) :: A
real*8, dimension(:), intent(in) :: v
real*8, dimension(size(v))  :: multmv 
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

logical function check_solve(A, x, b, eps)
real*8, dimension(:,:), intent(in) :: A
real*8, dimension(:), intent(in) :: x, b
real*8, intent(in), optional :: eps
real*8 :: e = 1.0E-5
real*8 :: summ
integer :: i,j,k,n
logical :: flag
check_solve = .true.
n = size(A, dim = 1)
if (presented(eps)) e = eps

if (check_solve) then
    !$OMP PARALLEL DO PRIVATE(i,j,summ) SHARED(A,x,b,check_solve,e,n)
    do i = 1, n
        summ = b(i)
        do j = 1, n
            summ = summ - A(i,j) * x(j)
        enddo
        if (abs(summ) > e) then
            !$OMP CRITICAL
             check_solve = .false.
            !$OMP END CRITICAL
        endif
    enddo
    !$OMP END PARALLEL DO
endif
end function check_solve

end module linalg
