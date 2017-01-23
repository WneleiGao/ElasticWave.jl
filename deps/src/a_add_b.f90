subroutine ca_add_b_rr(nthreads, n, a, b, c)
!DIR$ ATTRIBUTES DLLEXPORT :: ca_add_b_rr
!DIR$ ATTRIBUTES ALIAS: 'ca_add_b_rr_':: ca_add_b_rr
implicit none

integer(kind=8), intent(in) :: nthreads, n
real(kind=8)   , intent(in) :: a(n), b(n)
real(kind=8)   , intent(inout) :: c(n)

integer batchsize, bi, i, L
batchsize = 2046

!$OMP parallel do default(none), num_threads(nthreads),  &
!$OMP&         private(i, bi, L),                        &
!$OMP&         shared(n, batchsize, c, a, b)
do bi = 1, n, batchsize
   L = min(bi+batchsize-1, n)
   do i = bi, L
      c(i) = a(i) + b(i)
   end do
end do
!$OMP end parallel do

return
end subroutine ca_add_b_rr
