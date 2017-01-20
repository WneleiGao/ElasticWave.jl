using PyPlot, ElasticWave
nz = 177; nx = 203; ext= 30;   iflag = 1;
dx = 5. ; dz = 5. ; dt = 5e-4; tmax = 0.5; f0=30.;
vp =4000.*ones(nz,nx);  vs=2000.*ones(nz,nx); rho=2.5*ones(nz,nx);
fidMtx = CreateFidMtx(nz, nx, ext, iflag, vp, vs, rho, dz, dx, dt, f0);

root = homedir();
path_pv = join([root "/Desktop/pv.bin"])
isz = 2; isx = 102; ot = 0.0;
flags = vec([false false true true false]);
src = InitSource(isz, isx, nz, nx, ext, iflag, f0, ot, dt, flags);
MultiStepForward(path_pv, src, fidMtx, dz, dx)

path_wfd = join([root "/Desktop/wfd.bin"])
MultiStepForward(path_wfd, src, fidMtx)

wfd = ReadWfd(path_wfd, 500); (vxx, vxz, vzx, vzz) = ReadPv(path_pv, 500);
tmp = reshape(wfd.Vx, 177+30, 203+60); imshow(tmp[1:nz, ext+1:nx+ext]);
tmp = reshape(wfd.Vz, 177+30, 203+60); figure(); imshow(tmp[1:nz, ext+1:nx+ext]);
figure(); imshow(vzz+vxx); figure(); imshow(vzx-vxz);

dm = zeros(nz, nx); du = zeros(nz, nx); dm[88, 1:end] = 500.; du[88, 1:end] = 500.;
irx = collect(1:5:nx); irz = 2*ones(Int64, length(irx));
shotv = MultiStepForward(irz, irx, path_pv,  dm, du, fidMtx; tmax = 0.75)
imshow(shotv.Vx, aspect=0.04); figure(); imshow(shotv.Vz, aspect=0.04)

(dm, du) = MultiStepAdjoint(shotv, fidMtx, path_pv);
imshow(dm); figure(); imshow(du);
