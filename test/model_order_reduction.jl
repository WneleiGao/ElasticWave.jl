using PyPlot, ElasticWave
nz = 179; nx = 237; ext= 30;   iflag = 2;
dx = 20. ; dz = 20. ; dt = 2e-3; tmax = 2.0; f0=12.;
vp=4000.*ones(nz,nx);  vs=2000.*ones(nz,nx); rho=2.5e-4*ones(nz,nx);
fidMtx = CreateFidMtx(nz, nx, ext, iflag, vp, vs, rho, dz, dx, dt, f0);
# multi-sources case
isz = 1; isx = 118; ot = 0.0;
flags = vec([false false true true false]);
src = InitSource(isz, isx, nz, nx, ext, iflag, f0, ot, dt, flags);
# =========finite difference modeling==========
root = homedir();
path = join([root, "/Desktop/spt.bin"])
MultiStepForward(path, src, fidMtx, st=3*dt, otype="spt", tmax=tmax)

path1 = join([root, "/Desktop/base.bin"])
lambda = QRbase(path1, path, rk=200, style="mem")
Q = formStateMatrix(path1)
path_Fr = join([root, "/Desktop/Fr.bin"])
lowRankFd(path_Fr, fidMtx, path1, rk=200)
Q = formBaseMatrix(path1, rk)
path_FQ = join([root, "/Desktop/FQ"])
(nz, nx, ext, iflag, dt, nt) = InfoSnapShots(path_FQ)
if iflag == 1
   Nz = nz + 2*ext
elseif iflag == 2
   Nz = nz +   ext
end
Nx = nx + 2*ext
fid = open(path_FQ, "r"); seek(fid, sizeof(Float32)*5)
FQ = reshape(read(fid, Float32, Nz*Nx*10*nt), Nz*Nx*10, nt)
Fr = Q'*FQ
Fr1 = ReadRDFD(path_Fr)
vecnorm(Fr1-Fr)

# ===============test source projection(checked) =================
root = homedir()
path = join([root, "/Desktop/base.bin"])
spj = SrcProj(src, path)
Q = formBaseMatrix(path);
Qt = Q';
Nz = nz+ext; Nx=nx+2*ext; iz=isz; ix=isx+ext;
sc = zeros(size(Qt,1), src.nt);
for it = 1 : src.nt
    ind = (ix-1)*Nz+iz + Nz*Nx*4
    sc[:,it] = sc[:,it] + Qt[:,ind]*src.Txx[it]/2
    ind = (ix-1)*Nz+iz + Nz*Nx*5
    sc[:,it] = sc[:,it] + Qt[:,ind]*src.Txx[it]/2
    ind = (ix-1)*Nz+iz + Nz*Nx*6
    sc[:,it] = sc[:,it] + Qt[:,ind]*src.Tzz[it]/2
    ind = (ix-1)*Nz+iz + Nz*Nx*7
    sc[:,it] = sc[:,it] + Qt[:,ind]*src.Tzz[it]/2
end
vecnorm(spj.sc-sc)

# ========= test more MOR =======================
using ElasticWave
nz = 179; nx = 237; ext= 30;   iflag = 1;
dx = 20. ; dz = 20. ; dt = 2e-3; tmax = 2.0; f0=12.;
vp=4000.*ones(nz,nx);  vs=2000.*ones(nz,nx); rho=2.5e-4*ones(nz,nx);
fidMtx = CreateFidMtx(nz, nx, ext, iflag, vp, vs, rho, dz, dx, dt, f0);
# multi-sources case
isz = 1; isx = 118; ot = 0.0;
flags = vec([true false false false false]);
src = InitSource(isz, isx, nz, nx, ext, iflag, f0, ot, dt, flags);
# =========finite difference modeling==========
root = homedir();
path = join([root, "/Desktop/spt.bin"])
MultiStepForward(path, src, fidMtx, st=dt, otype="spt", tmax=tmax)

rk = 200
pathe = join([root, "/Desktop/base.bin"])
lambda = QRbase(pathe, path, rk=rk, style="mem")

pathFr = join([root, "/Desktop/Fr.bin"])
lowRankFd(pathFr, fidMtx, pathe;)

spj = SrcProj(src, pathe, rk=rk)
pathc = join([root, "/Desktop/wc.bin"])
Fr = ReadRDFD(pathFr);

MultiStepForward_MOR(pathc, pathFr, spj, tmax=tmax)
path1 = join([root, "/Desktop/spt1.bin"])
ModSpts2Spts(path1, pathe, pathc, rk=rk)

path1 = join([root, "/Desktop/spt1.bin"])
MultiStepForward(path1, src, fidMtx, st=dt, otype="spt", tmax=tmax)

Nz = nz+ext; Nx = nx+2*ext; len = Nz*Nx;
Vx = tmp[1:len] + tmp[len+1:2*len]
Vx = reshape(Vx, Nz, Nx)





spt = ReadSnapShot(pathe, 50);
vx = reshape(spt.Vxx+spt.Vxz, 239, 297);
SeisPlot(vx)

spt = ReadSnapShot(pathe, 120);
txx = reshape(spt.Txxx+spt.Txxz, 239, 297);
SeisPlot(txx)












# path_spt2vx = join([root, "/Desktop/spt2vx.bin"])
# ExtractCpt(path_spt2vx, path, itype="spt")
# path_wfd2vx = join([root, "/Desktop/wfd2vx.bin"])
# ExtractCpt(path_wfd2vx, path1, itype="wfd")
#
# fid  = open(path_spt2vx, "r")
# fid1 = open(path_wfd2vx, "r")
# tmp = 0.0; tmp1 = 0.0
# nz=read(fid, Int32); nx=read(fid, Int32); ext=read(fid, Int32); iflag=read(fid, Int32); dt=read(fid, Float32);
# sptSize = sizeof(Float32)*nz*nx
# nt = floor(Int64, (filesize(fid)-sizeof(Float32)*5)/sptSize)
#     it = 450
#     pos = (it-1)*sptSize +  sizeof(Float32)*5
#     seek(fid,  pos)
#     seek(fid1, pos)
#     vx = read(fid, Float32, nz*nx)
#     vx1 = read(fid1, Float32, nz*nx)
#     vx = reshape(vx, Int64(nz), Int64(nx))
#     figure()
#     vx1 = reshape(vx1, Int64(nz), Int64(nx))
