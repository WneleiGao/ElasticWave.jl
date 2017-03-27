subroutine a_mul_b_prr(nthreads, m, n, A, jA, iA, x, y, yt)

!DIR$ ATTRIBUTES DLLEXPORT :: a_mul_b_prr
!DIR$ ATTRIBUTES ALIAS: 'a_mul_b_prr_':: a_mul_b_prr

use omp_lib
implicit none

integer(kind=8),intent(in)   :: m, n, nthreads  ! # of rows in A
real(kind=8)   ,intent(in)   :: A(*), jA(*), iA(n+1), x(n)
real(kind=8)   ,intent(inout):: y(m), yt(m*nthreads)

integer i, j1, j2, j, jaj, mythread, mm, jm
real(kind=8) xi

!$OMP parallel default(none), num_threads(nthreads), &
!$OMP&         private(mythread, i, xi, j1,j2, j, jaj, mm, jm),  &
!$OMP&         shared(m,n, A,jA,iA, x,y, yt, nthreads)

mythread = OMP_GET_THREAD_NUM()
mm = m * mythread
yt(mm+1 : mm+m) = 0.d0

!$OMP do
do i = 1, n
   xi = x(i)
   j1 = iA(i)
   j2 = iA(i+1) - 1
   do j = j1, j2
      jaj = jA(j)
      yt(mm + jaj) = yt(mm + jaj) + xi*A(j)
    end do  ! j
end do  ! i
!$OMP end do

do j = 0, nthreads-1
   jm = j*m
   !$OMP do
   do i = 1, m
      y(i) = y(i) + yt(i+jm)
   end do
   !$OMP end do
end do ! j

!$OMP end parallel

return
end subroutine a_mul_b_prr

! =======================real matrix * real vector================================
subroutine a_mul_b_rr(m, n, A, jA, iA, x, y)

!DIR$ ATTRIBUTES DLLEXPORT :: a_mul_b_rr
!DIR$ ATTRIBUTES ALIAS: 'a_mul_b_rr_':: a_mul_b_rr

implicit none

integer(kind=8),intent(in)   :: m  ! # of rows in A
integer(kind=8),intent(in)   :: n  ! # of columns in A
real(kind=8)   ,intent(in)   :: A(*)
integer(kind=8),intent(in)   :: jA(*), iA(n+1)
real(kind=8)   ,intent(inout):: y(m)
real(kind=8)   ,intent(in)   :: x(n)

integer i, j1, j2, j, jaj

include "a_mul_b.fi"

return
end subroutine a_mul_b_rr

! =======================real matrix * complex vector================================
subroutine a_mul_b_rc(m, n, A, jA, iA, x, y)

!DIR$ ATTRIBUTES DLLEXPORT :: a_mul_b_rc
!DIR$ ATTRIBUTES ALIAS: 'a_mul_b_rc_':: a_mul_b_rc

implicit none

integer(kind=8),intent(in)   :: m  ! # of rows in A
integer(kind=8),intent(in)   :: n  ! # of columns in A
real(kind=8)   ,intent(in)   :: A(*)
integer(kind=8),intent(in)   :: jA(*), iA(n+1)
complex(kind=8),intent(inout):: y(m)
complex(kind=8),intent(in)   :: x(n)

integer i, j1, j2, j, jaj

include "a_mul_b.fi"

return
end subroutine a_mul_b_rc

!========================complex matrix * complex vector============================================
subroutine a_mul_b_cc(m, n, A, jA, iA, x, y)

!DIR$ ATTRIBUTES DLLEXPORT :: a_mul_b_cc
!DIR$ ATTRIBUTES ALIAS: 'a_mul_b_cc_':: a_mul_b_cc

implicit none

integer(kind=8),intent(in)   :: m  ! # of rows in A
integer(kind=8),intent(in)   :: n  ! # of columns in A
complex(kind=8),intent(in)   :: A(*)
integer(kind=8),intent(in)   :: jA(*), iA(n+1)
complex(kind=8),intent(inout):: y(m)
complex(kind=8),intent(in)   :: x(n)

integer i, j1, j2, j, jaj

include "a_mul_b.fi"

return
end subroutine a_mul_b_cc
