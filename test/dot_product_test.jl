using PyPlot, ElasticWave

nz = 200; nx = 200; ext= 30;   iflag = 1;
dx = 5. ; dz = 5. ; dt = 5e-4; tmax = 0.5; f0=30.;
vp=4000.*ones(nz,nx);  vs=2000.*ones(nz,nx); rho=2.5*ones(nz,nx);
fidMtx = CreateFidMtx(nz, nx, ext, iflag, vp, vs, rho, dz, dx, dt, f0);

# multi-sources case
isz = 2; isx = 100; ot = 0.0;
flags = vec([true false false false false]);
src = InitSource(isz, isx, nz, nx, ext, iflag, f0, ot, dt, flags);

# =========finite difference modeling==========
path_pv = "/Users/wenlei/Desktop/pv.bin"
MultiStepForward(path_pv, src, fidMtx, dz, dx);
# (vxx, vxz, vzx, vzz) = ReadPv(path, 620); imshow(vxx+vzz); figure(); imshow(vxz-vzx);

path = "/Users/wenlei/Desktop/test.bin"
irx = collect(1:5:nx); irz = 2*ones(Int64, length(irx));
dm  = zeros(nz, nx); dm[100,100] = 0.5;
du  = zeros(nz, nx); du[100,100] = 0.5;
shotv = MultiStepForward(irz, irx, path, dm, du, fidMtx, tmax=0.75);
path_shotv = "/Users/wenlei/Desktop/shotv.bin"
WriteShotV(path_shotv)


using PyPlot, ElasticWave
nz = 200; nx = 200; ext= 30;   iflag = 1;
dx = 5. ; dz = 5. ; dt = 5e-4; tmax = 0.5; f0=30.;
vp=4000.*ones(nz,nx);  vs=2000.*ones(nz,nx); rho=2.5*ones(nz,nx);
fidMtx = CreateFidMtx(nz, nx, ext, iflag, vp, vs, rho, dz, dx, dt, f0);
path = "/Users/wenlei/Desktop/test.bin"
path_shotv = "/Users/wenlei/Desktop/shotv.bin"
shotv = ReadShotV(path_shotv)
path_wfd = "/Users/wenlei/Desktop/wfd_recon.bin"
MultiStepAdjoint(path_wfdrecon, shotv, fidMtx)

# wfd = ReadWfd(path_wfdrecon, 200); imshow(reshape(wfd.Vx, 260, 260))
# figure(); wfd = ReadWfd(path_wfdrecon, 300); imshow(reshape(wfd.Vx, 260, 260))
(dm1, du1) = MultiStepAdjoint(shotv, fidMtx, path)
