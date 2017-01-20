# test one step forward and adjoint
using PyPlot, ElasticWave

nz = 233; nx = 277; ext= 30;   iflag = 2;
dx = 5. ; dz = 5. ; dt = 5e-4; tmax = 0.5; f0=30.;
vp =4000.*ones(nz,nx);  vs=2000.*ones(nz,nx); rho=2.5*ones(nz,nx);
fidMtx = CreateFidMtx(nz, nx, ext, iflag, vp, vs, rho, dz, dx, dt, f0);

isz = 2; isx = 138; ot = 0.0;
flags = vec([true false false false false]);
src = InitSource(isz, isx, nz, nx, ext, iflag, f0, ot, dt, flags);

spt1 = InitSnapShot(nz, nx, ext, iflag, dt, 1);
spt2 = InitSnapShot(nz, nx, ext, iflag, dt, 2);

spt1.Vxx = randn(length(spt1.Vxx));
spt1.Vxz = randn(length(spt1.Vxz));
spt1.Vzx = randn(length(spt1.Vzx));
spt1.Vzz = randn(length(spt1.Vzz));
spt1.Txxx = randn(length(spt1.Txxx));
spt1.Txxz = randn(length(spt1.Txxz));
spt1.Tzzx = randn(length(spt1.Tzzx));
spt1.Tzzz = randn(length(spt1.Tzzz));
spt1.Txzx = randn(length(spt1.Txzx));
spt1.Txzz = randn(length(spt1.Txzz));
OneStepForward!(spt2, spt1, fidMtx);

spt4 = InitSnapShot(nz, nx, ext, iflag, dt, 4);
spt3 = InitSnapShot(nz, nx, ext, iflag, dt, 3);
spt4.Vxx = randn(length(spt1.Vxx));
spt4.Vxz = randn(length(spt1.Vxz));
spt4.Vzx = randn(length(spt1.Vzx));
spt4.Vzz = randn(length(spt1.Vzz));
spt4.Txxx = randn(length(spt4.Txxx));
spt4.Txxz = randn(length(spt4.Txxz));
spt4.Tzzx = randn(length(spt4.Tzzx));
spt4.Tzzz = randn(length(spt4.Tzzz));
spt4.Txzx = randn(length(spt4.Txzx));
spt4.Txzz = randn(length(spt4.Txzz));
OneStepAdjoint!(spt3, spt4, fidMtx);

tmp1 = IptSpts(spt1, spt3)
tmp2 = IptSpts(spt4, spt2)

(vxx, vxz, vzx, vzz) = ReadPv(path_pv, 500)
imshow(vxx+vzz); figure();  imshow(vzx-vxz);
