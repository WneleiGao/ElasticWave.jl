
using PyPlot, ElasticWave

nz = 123; nx = 321; ext= 20;   iflag = 1;
dx = 5. ; dz = 5. ; dt = 5e-4; tmax = 0.5; f0=30.;
vp =4000.*ones(nz,nx);  vs=2000.*ones(nz,nx); rho=2.5*ones(nz,nx)*1e-4;
fidMtx = CreateFidMtx(nz, nx, ext, iflag, vp, vs, rho, dz, dx, dt, f0);

# test one step forward and adjoint
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

tmp1 = IptSpts(spt1, spt5)
tmp2 = IptSpts(spt2, spt4)
# =========================pass one-step dot product test=======================

# multi-step dot product test
isz = []; isx = []; ot = [];
flags = [false, false, true, true, false]
for iz = 1 : nz
    for ix = 1 : nx
        if rand() < 0.01
           push!(isz, iz)
           push!(isx, ix)
           push!(ot, floor(Int64, 0.5*rand()/dt)*dt)
        end
    end
end
isz = convert(Array{Int64,1}, isz);
isx = convert(Array{Int64,1}, isx);
ot  = convert(Array{Float64,1}, ot);
srcs = InitMultiSources(isz, isx, nz, nx, ext, iflag, f0, ot, dt, flags);
root = homedir(); path = joinpath(root, "Desktop/spt.bin");
srcs2spt(path, srcs);
irx = collect(1:nx); irz = ones(Int64, length(irx));
shotv = MultiStepForward(irz, irx, path, fidMtx, tmax=1.0);

nt = shotv.nt
shotv_adj = InitShotV(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt);
wlet = Ricker(f0, dt); nw = length(wlet);
for ir = 1 : length(irz)
    if rand() < 0.5
       tmp = -1. * wlet
    else
       tmp =  wlet
    end
    inds = floor(Int64, rand()*(nt-nw))
    shotv_adj.Vx[inds:inds+nw-1 ,ir] = tmp[:]
    shotv_adj.Vz[inds:inds+nw-1 ,ir] = tmp[:]
end
path_adj = joinpath(root, "Desktop/adj.bin")
MultiStepAdjoint_spt(path_adj, shotv_adj, fidMtx)

tmp = dot(vec(shotv.Vx), vec(shotv_adj.Vx))+dot(vec(shotv.Vz), vec(shotv_adj.Vz))
tmp1 = 0.
(nz, nx, ext, iflag, dt, nt) = InfoSnapShots(path)
for it = 1 : nt
    spt1 = ReadSnapShot(path, it)
    spt2 = ReadSnapShot(path_adj, it)
    tmp1 = tmp1 + IptSpts(spt1, spt2)
end
# =========================pass multi-step dot product test=====================

# test born approximation
isz = 1; isx = floor(Int64, nx/2); ot = 0.0;
flags = vec([false false true true false]);
src = InitSource(isz, isx, nz, nx, ext, iflag, f0, ot, dt, flags);
root = homedir(); path = joinpath(root, "Desktop/pv.bin")
MultiStepForward(path, src, fidMtx, dz, dx, tmax=1.0)
if iflag == 1
   Nz = nz + 2*ext
elseif iflag == 2
   Nz = nz +   ext
end
Nx = nx + 2*ext
dm = vec(100*randn(Nz, Nx)); du = vec(100*randn(Nz, Nx));
irx = collect(1:nx); irz = ones(Int64, length(irx));
shotv = MultiStepForward(irz, irx, path, dm, du, fidMtx, tmax=1.0)

nt = shotv.nt
adj = InitShotV(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt);
wlet = Ricker(f0, dt); nw = length(wlet);
for ir = 1 : length(irz)
    if rand() < 0.5
       tmp = -1. * wlet
    else
       tmp =  wlet
    end
    inds = floor(Int64, rand()*(nt-nw))
    adj.Vx[inds:inds+nw-1 ,ir] = tmp[:]
    adj.Vz[inds:inds+nw-1 ,ir] = tmp[:]
end
(dm1, du1) = MultiStepAdjoint(adj, fidMtx, path)
tmp = dot(vec(shotv.Vx), vec(adj.Vx))+dot(vec(shotv.Vz), vec(adj.Vz))
tmp1 = dot(dm, dm1)+dot(du, du1)
# ==============born approximation=======================

# test one-shot imaging
isz = 1; isx = floor(Int64, nx/2); ot = 0.0;
flags = vec([false false true true false]);
src = InitSource(isz, isx, nz, nx, ext, iflag, f0, ot, dt, flags);
root = homedir(); path = joinpath(root, "Desktop/pv.bin")
MultiStepForward(path, src, fidMtx, dz, dx, tmax=1.0)
if iflag == 1
   Nz = nz + 2*ext
elseif iflag == 2
   Nz = nz +   ext
end
Nx = nx + 2*ext
dm = zeros(Nz, Nx); dm[80,:]=500; dm = vec(dm);
du = zeros(Nz, Nx); du[80,:]=500; du = vec(du);
irx = collect(1:nx); irz = ones(Int64, length(irx));
shotv = MultiStepForward(irz, irx, path, dm, du, fidMtx, tmax=1.0)

(dm1, du1) = MultiStepAdjoint(shotv, fidMtx, path);
SeisPlot(reshape(dm1, Nz, Nx));
SeisPlot(reshape(du1, Nz, Nx));

(vxx, vxz, vzx, vzz) = ReadPv(path, 300);
p = reshape(vxx+vzz, Nz, Nx); SeisPlot(p);
s = reshape(vxz-vzx, Nz, Nx); SeisPlot(s);
