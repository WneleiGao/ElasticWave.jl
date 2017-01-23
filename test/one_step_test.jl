# test one step forward and adjoint
using PyPlot, ElasticWave

nz = 300; nx = 300; ext= 30;   iflag = 2;
dx = 5. ; dz = 5. ; dt = 5e-4; tmax = 0.5; f0=30.;
vp =4000.*ones(nz,nx);  vs=2000.*ones(nz,nx); rho=2.5*ones(nz,nx);
fidMtx = CreateFidMtx(nz, nx, ext, iflag, vp, vs, rho, dz, dx, dt, f0);

isz = 2; isx = 138; ot = 0.0;
flags = vec([true false false false false]);
src = InitSource(isz, isx, nz, nx, ext, iflag, f0, ot, dt, flags);

spt1 = InitSnapShot(nz, nx, ext, iflag, dt, 1);
spt2 = InitSnapShot(nz, nx, ext, iflag, dt, 2);
spt3 = InitSnapShot(nz, nx, ext, iflag, dt, 2);

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
tmp  = zeros(length(spt1.Txzz));
tmp1 = zeros(length(spt1.Txzz));
@time OneStepForward!(spt2, spt1, fidMtx, tmp, tmp1);
@time OneStepForward!(spt3, spt1, fidMtx);
d = dif2spt(spt2, spt3)

spt4 = InitSnapShot(nz, nx, ext, iflag, dt, 2);
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
spt5 = InitSnapShot(nz, nx, ext, iflag, dt, 1);
spt6 = InitSnapShot(nz, nx, ext, iflag, dt, 2);

@time OneStepAdjoint!(spt5, spt4, fidMtx);
@time OneStepAdjoint!(spt6, spt4, fidMtx, tmp, tmp1);
d = dif2spt(spt5, spt6)

function IptSpts_test(spt1::SnapShot, spt2::SnapShot)
    d = 0.0
    d = d + dot(spt1.Vxx, spt2.Vxx)
    d = d + dot(spt1.Vxz, spt2.Vxz)
    d = d + dot(spt1.Vzx, spt2.Vzx)
    d = d + dot(spt1.Vzz, spt2.Vzz)
    d = d + dot(spt1.Txxx, spt2.Txxx)
    d = d + dot(spt1.Txxz, spt2.Txxz)
    d = d + dot(spt1.Tzzx, spt2.Tzzx)
    d = d + dot(spt1.Tzzz, spt2.Tzzz)
    d = d + dot(spt1.Txzx, spt2.Txzx)
    d = d + dot(spt1.Txzz, spt2.Txzz)
    return d
end

tmp1 = IptSpts(spt1, spt5)
tmp2 = IptSpts(spt2, spt4)
