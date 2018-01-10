! update for sigma_zz, an sigma_xx
subroutine forward(nz, nx, nstep, &
                   deltaZ, deltaX, dt, &
                   lambda, mu, rho, &
                   hLambdaZ, hLambda2MuZ, hMuZ, hMuX, hRhoZX, &
                   sigmazz, sigmaxx, sigmazx, vz, vx, &
                   vdvzdz, vdvxdx, vdvxdz, vdvzdx, vdsigmazzdz, vdsigmazxdx, vdsigmaxxdx, vdsigmazxdz, &
                   mdvzdz, mdvxdx, mdvxdz, mdvzdx, mdsigmazzdz, mdsigmazxdx, mdsigmaxxdx, mdsigmazxdz, &
                   az, bz, kz, haz, hbz, hkz, &
                   ax, bx, kx, hax, hbx, hkx   )
    implicit none

    integer(kind=8)                         :: nz, nx, nstep
    real(kind=8)                            :: deltaZ, deltaX, dt
    real(kind=8), dimension(nz, nx)         :: lambda, mu, rho
    real(kind=8)                            :: hLambdaZ, hLambda2MuZ, hMuz, hMuX, hRhoZX
    real(kind=8), dimension(0:nz+1, 0:nx+1) :: sigmazz, sigmaxx, sigmazx, vz, vx
    real(kind=8)                            :: vdvzdz, vdvxdx, vdvxdz, vdvzdx, vdsigmazzdz, vdsigmazxdx, vdsigmaxxdx, vdsigmazxdz
    real(kind=8), dimension(nz, nx)         :: mdvzdz, mdvxdx, mdvxdz, mdvzdx, mdsigmazzdz, mdsigmazxdx, mdsigmaxxdx, mdsigmazxdz

    real(kind=8), dimension(nz)             :: az, bz, kz, haz, hbz, hkz
    real(kind=8), dimension(nx)             :: ax, bx, kx, hax, hbx, hkx

    integer(kind=8)                         :: i, j, it
    ! source related terms
    real(kind=8) :: f0, pi, factor, angle, t0, t, a, source, forceZ, forceX

    do it = 1 , nstep

      ! update sigmazz and sigmaxx
      do i = 2 , nx
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

      ! update sigmazx
      do j = 1 , nx-1
         do i = 2 , nz
            hMuX = 0.5d0 * (mu(i,j+1) + mu(i,j))

            vdvxdz = (27.d0*vx(i,j  ) - 27.d0*vx(i-1,j) - vx(i+1,j  ) + vx(i-2,j  )) / (24.d0*deltaZ)
            vdvzdx = (27.d0*vz(i,j+1) - 27.d0*vz(i  ,j) - vz(i  ,j+2) + vz(i  ,j-1)) / (24.d0*deltaX)
            mdvxdz(i,j) =  bz(i) * mdvxdz(i,j) +  az(i) * vdvxdz
            mdvzdx(i,j) = hbx(j) * mdvzdx(i,j) + hax(j) * vdvzdx

            vdvxdz = vdvxdz /  kz(i) + mdvxdz(i,j)
            vdvzdx = vdvzdx / hkx(j) + mdvzdx(i,j)  !!double check the damping parameter!!

            sigmazx(i,j) = sigmazx(i,j) + hMuX * (vdvxdz + vdvzdx) * dt
         enddo
      enddo

      ! update vz
      do j = 2 , nx
         do i = 2 , nz
            vdsigmazzdz = (27.d0*sigmazz(i,j) - 27.d0*sigmazz(i-1,j) - sigmazz(i+1,j) + sigmazz(i-2,j)) / (24.d0*deltaZ)
            vdsigmazxdx = (27.d0*sigmazx(i,j) - 27.d0*sigmazx(i,j-1) - sigmazx(i,j+1) + sigmazx(i,j-2)) / (24.d0*deltaX)

            mdsigmazzdz(i,j) = bz(i) * mdsigmazzdz(i,j) + az(i) * vdsigmazzdz
            mdsigmazxdx(i,j) = bx(j) * mdsigmazxdx(i,j) + ax(j) * vdsigmazxdx

            vdsigmazzdz = vdsigmazzdz / kz(i) + mdsigmazzdz(i,j)
            vdsigmazxdx = vdsigmazxdx / kx(j) + mdsigmazxdx(i,j)

            vz(i,j) = vz(i,j) + (vdsigmazzdz + vdsigmazxdx) * dt / rho(i,j)
         enddo
      enddo

      ! update vx
      do j = 1 , nx-1
         do i = 1 , nz-1
            hRhoZX = 0.25d0 * (rho(i,j) + rho(i+1,j) + rho(i,j+1) + rho(i+1,j+1))

            vdsigmazxdz = (27.d0*sigmazx(i+1,j) - 27.d0*sigmazx(i,j) - sigmazx(i+2,j) + sigmazx(i-1,j)) / (24.d0*deltaZ)
            vdsigmaxxdx = (27.d0*sigmaxx(i,j+1) - 27.d0*sigmaxx(i,j) - sigmaxx(i,j+2) + sigmazx(i,j-1)) / (24.d0*deltaX)

            mdsigmazxdz(i,j) = hbz(i) * mdsigmazxdz(i,j) + haz(i) * vdsigmazxdz
            mdsigmaxxdx(i,j) = hbx(j) * mdsigmaxxdx(i,j) + hax(j) * vdsigmaxxdx

            vdsigmazxdz = vdsigmazxdz / hkz(i) + mdsigmazxdz(i,j)
            vdsigmaxxdx = vdsigmaxxdx / hkx(j) + mdsigmaxxdx(i,j)
         enddo
      enddo

      ! add source to wave field
      f0 = 7.0; pi = 3.1415926d0; factor = 1.d7; angle = 2.3561d0
      t0 = 1.2d0 / f0
      a = pi * pi * f0 * f0
      t = (it-1)  * dt
      source = - factor * 2.d0 * a * (t-t0) * exp(-a*(t-t0)**2)
      forceZ = sin(angle) * source
      forceX = cos(angle) * source

      i = 20; j = 320;
      hRhoZX = 0.25d0 * (rho(i,j) + rho(i+1,j) + rho(i,j+1) + rho(i+1,j+1))

      vz(i,j) = vz(i,j) + forceZ * dt / rho(i,j)
      vx(i,j) = vx(i,j) + forceX * dt / hRhoZX
    enddo

end subroutine forward



subroutine assign(A, m, n)
    implicit none
    integer(kind=8), intent(in) :: m, n
    integer(kind=8), intent(inout), dimension(0:m+1, 0:n+1) :: A

    integer(kind=8) :: i, j, c=1
    real(kind=8) :: v

    do j = 1, n
       do i = 1, m
          A(i,j) = c
          c = c + 1
       end do
    end do

    A(1,1) = 1.d7

end subroutine assign
