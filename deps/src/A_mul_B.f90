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
