! update for sigma_zz, an sigma_xx
subroutine OneTimeStepForward(nz, nx, &
                              sigmazz, sigmaxx, sigmaxz, vz, vx, &
                              hLambda2Muz, vdvzdz, mdvzdz, &
                              hMux       , vdvxdx, mdvxdx   )

integer(kind=8) :: nz, nx
real(kind=8) :: hLambda2Muz, vdvzdz, mdvzdz, hMux, vdvxdx, mdvxdx
real(kind=8), dimension(0:nz+1, 0:nx+1) :: sigmazz, sigmaxx, sigmaxz, vz, vx

integer(kind=4) :: i, j

do i = 2 : nx
   do j = 1, nz

   enddo
enddo

end subroutine OneTimeStepForward
