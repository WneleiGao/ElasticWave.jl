! update for sigma_zz, an sigma_xx
subroutine OneTimeStepForward(nz, nx, &
                              sigmazz, sigmaxx, sigmaxz, vz, vx, &
                              hLambda2Muz, vdvzdz, mdvzdz, &
                              hMux       , vdvxdx, mdvxdx   )

    integer(kind=8) :: nz, nx
    real(kind=8) :: hLambdaZ, hLambdaZ, hLambda2MuZ, vdvzdz, mdvzdz, hMuX, vdvxdx, mdvxdx
    real(kind=8), dimension(0:nz+1, 0:nx+1) :: sigmazz, sigmaxx, sigmaxz, vz, vx

    integer(kind=4) :: i, j

    do i = 2 : nx
       do j = 1, nz
          hLambdaZ    = 0.5d0 * (lambda(i,j) + lambda(i+1,j))
          hMuZ        = 0.5d0 * (    mu(i,j) +     mu(i+1,j))
          hLambda2MuZ = hLambdaZ + 2.d0 * hmuz
          ! partial derivative
          vdvzdz = (27.d0*vz(i+1,j)-27.d0*vz(i,j  ) - vz(i+2,j) + vz(i-1,j)) / (24.d0*deltaZ)
          vdvxdx = (27.d0*vx(i  ,j)-27.d0*vx(i,j-1) - vx(i,j+1) + vx(i,j-2)) / (24.d0*deltaX)
          ! update memory variable
          mdvzdz(i,j) = hbz(i)*mdvzdz(i,j) + haz(i) * vdvzdz
          mdvxdx(i,j) =  bx(j)*mdvxdx(i,j) +  ax(j) * vdvxdx
          ! apply dampping
          vdvzdz = vdvzdz / hkz(i) + mdvzdz(i,j)
          vdvxdx = vdvxdx /  kx(j) + mdvxdx(i,j)
          ! update state variable
          sigmazz(i,j) = sigmazz(i,j) + (hLambda2MuZ * vdvzdz + hLambdaZ * vdvxdx) * dt
          sigmaxx(i,j) = sigmaxx(i,j) + (hLambda2MuZ * vdvxdx + hLambdaZ * vdvzdz) * dt
       enddo
    enddo

    do j = 1 , nx-1
       do i = 2 , nz
          hMuX = 0.5d0 * (mu(i,j+1) + mu(i,j))

          vdvxdz = (27.d0*vx(i,j  ) - 27.d0*vx(i-1,j) - vx(i+1,j  ) + vx(i-2,j  )) / (24.d0*deltaZ)
          vdvzdx = (27.d0*vz(i,j+1) - 27.d0*vz(i  ,j) - vz(i  ,j+2) + vz(i  ,j-1)) / (24.d0*deltaX)

          mdvxdz(i,j) =  bz(i) * mdvxdz(i,j) +  az(i) * vdvxdz(i,j)
          mdvzdx(i,j) = hbx(j) * mdvzdx(i,j) + hax(j) * vdvzdx(i,j)

          vdvxdz = vdvxdz /  kz(i) + mdvxdz(i,j)
          vdvzdx = vdvzdx / hkx(j) + mdvzdx(i,j)  !!double check the damping parameter!!
       enddo
    enddo

    do j = 2 , nx
       do i = 2 , nz
          vdsigmazzdz = (27.d0*sigmazz(i,j) - 27.d0*sigmazz(i-1,j) - sigmazz(i+1,j) + sigmazz(i-2,j)) / (24.d0*deltaZ)
          vdsigmazxdx = (27.d0*sigmazx(i,j) - 27.d0*sigmazx(i,j-1) - sigmazx(i,j+1) + sigmazx(i,j-2)) / (24.d0*deltaX)

          mdsigmazzdz(i,j) = bz(i) * mdsigmazzdz(i,j) + az(i) * vdsigmazzdz
          mdsigmazxdx(i,j) = bx(j) * mdsigmazxdx(i,j) + ax(j) * vdsigmazxdx

          vdsigmazzdz = vdsigmazzdz / kz(i) + mdsigmazzdz
          vdsigmazxdx = vdsigmazxdx / kx(j) + mdsigmazxdx

          
       enddo
    enddo


end subroutine OneTimeStepForward
