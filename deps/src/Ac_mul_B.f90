subroutine ac_mul_b_rr(nthreads, m, n, A, jA, iA, x, y)
!DIR$ ATTRIBUTES DLLEXPORT :: ac_mul_b_rr
!DIR$ ATTRIBUTES ALIAS: 'ac_mul_b_rr_':: ac_mul_b_rr

implicit none

integer(kind=8),intent(in)   :: nthreads
integer(kind=8),intent(in)   :: m  ! # of rows in A
integer(kind=8),intent(in)   :: n  ! # of columns in A

real(kind=8)   ,intent(in)   :: A(*)
integer(kind=8),intent(in)   :: jA(*), iA(n+1)
real(kind=8)   ,intent(in)   :: x(m)
real(kind=8)   ,intent(inout):: y(n)

integer i, j1, j2, j
real(kind=8) t

#include "ac_mul_b.fi"

return
end subroutine ac_mul_b_rr


!------------------------------------------------------------------------
subroutine ac_mul_b_rc(nthreads, m, n, A, jA, iA, x, y)
!DIR$ ATTRIBUTES DLLEXPORT :: ac_mul_b_rc
!DIR$ ATTRIBUTES ALIAS: 'ac_mul_b_rc_':: ac_mul_b_rc

implicit none

integer(kind=8),intent(in)   :: nthreads
integer(kind=8),intent(in)   :: m  ! # of rows in A
integer(kind=8),intent(in)   :: n  ! # of columns in A

real(kind=8)   ,intent(in)   :: A(*)
integer(kind=8),intent(in)   :: jA(*), iA(n+1)
complex(kind=8),intent(in)   :: x(m)
complex(kind=8),intent(inout):: y(n)

integer i, j1, j2, j
complex(kind=8) t

#include "ac_mul_b.fi"

return
end subroutine ac_mul_b_rc

!------------------------------------------------------------------------
subroutine ac_mul_b_cc(nthreads, m, n, A, jA, iA, x, y)
!DIR$ ATTRIBUTES DLLEXPORT :: ac_mul_b_cc
!DIR$ ATTRIBUTES ALIAS: 'ac_mul_b_cc_':: ac_mul_b_cc

implicit none

integer(kind=8),intent(in)   :: nthreads
integer(kind=8),intent(in)   :: m  ! # of rows in A
integer(kind=8),intent(in)   :: n  ! # of columns in A

complex(kind=8),intent(in)   :: A(*)
integer(kind=8),intent(in)   :: jA(*), iA(n+1)
complex(kind=8),intent(in)   :: x(m)
complex(kind=8),intent(inout):: y(n)

integer i, j1, j2, j
complex(kind=8) t

integer bi, L, rowBatchSize

rowBatchSize = 1536
y = 0.0
!$OMP parallel do default(none), num_threads(nthreads),  &
!$OMP&         private(i, t, j1, j2, j, bi, L),          &
!$OMP&         shared(rowBatchSize, n, A, jA, iA, x, y)
do bi = 1, n, rowBatchSize
	 L  = min(bi+rowBatchSize-1,n)
	 do i = bi,L
			j1 = iA(i)
			j2 = iA(i+1) - 1
			t  = 0.d0
			do  j = j1, j2
					t = t + conjg(A(j)) * x(jA(j))
		  end do !j
			y(i) = y(i) + t
	 end do !i
end do ! bi
!$OMP end parallel do

return
end subroutine ac_mul_b_cc
