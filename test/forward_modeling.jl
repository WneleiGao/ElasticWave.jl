using PyPlot, ElasticWave

nz = 300; nx = 300; # model size
ext= 20; # thickness of PML layers
iflag = 1; #when iflag=1, four side obsorbing boundary, iflag=2, free surface
dx = 5. ; dz = 5. ; dt = 4e-4; # grid size
tmax = 1.5; f0=20.; # recording time and dominant frequency of source wavelet
vp=2500.*ones(nz,nx); vs=2500./sqrt(3)*ones(nz,nx); rho=2.5*ones(nz,nx);
vp[151:end,:]=3500.;  vs[151:end,:]=3500./sqrt(3.); #physical model
SeisPlot(vp, xlabel="X (km)", ylabel="Z (km)", vmin=1500, vmax=4500, cmap="gray", dx=0.005, dy=0.005);
tight_layout(); savefig("/Users/wenyue/Desktop/v.pdf")
fidMtx = CreateFidMtx(nz, nx, ext, iflag, vp, vs, rho, dz, dx, dt, f0); #finite difference matrix

# single source
isz = 5; isx = 150; ot = 0.0; #source location and starting time
flags = vec([false false true true false]); # add source to normal stress, simulate explosive sources
src = InitSource(isz, isx, nz, nx, ext, iflag, f0, ot, dt, flags);
src.Txx = src.Txx * 1e5; src.Tzz = src.Tzz * 1e5;


# # =========finite difference modeling==========
# irx = collect(1:nx); irz = 1*ones(Int64, length(irx)); #location of receivers
# shotv = MultiStepForward(irz, irx, src, fidMtx, tmax=tmax) #generate one shot
root = homedir(); path = join([root "/Desktop/wfd.bin"]);
MultiStepForward(path, src, fidMtx, tmax=tmax);

root = homedir(); path = join([root "/Desktop/decomposed.bin"]);
MultiStepForward(path, src, fidMtx, dz, dx; tmax=tmax)

root = homedir(); path = join([root "/Desktop/wfd.bin"]);
MultiStepForward(path, src, fidMtx, tmax=tmax);

(nz, nx, ext, iflag, dt, nt) = InfoPv(path);
pp = zeros(nt); ps = zeros(nt);
rz = irz[1] + ext; rx = irx[1] + ext; idx = (rx-1)*(nz+2*ext) + rz + ext;
for i = 1151 : nt
    (vxx, vxz, vzx, vzz) = ReadPv(path, i)
    pp[i] = vxx[idx] + vzz[idx]
    ps[i] = vxz[idx] - vzx[idx]
end

path1 = join([root "/Desktop/moive"]);
waveAnim(path1, path, cpt="Vz", cmap="gray", wbox=6, hbox=6, dt1=0.003, pclip=1.0, vlim=0.1)

irz = [1]; irx = [250];
shotv = binWfd2shotv(irz, irx, path)

obs = hcat(shotv.Vz, shotv.Vx); obs[1:1150,:] = 0.0;
spt = hcat(pp, ps);
SeisPlot(obs, style="wiggles", dy=0.0004, xcur=0.8, ylabel="Time (s)", xticks=[]);
tight_layout(); savefig("/Users/wenyue/Desktop/obs.pdf");

SeisPlot(spt, style="wiggles", dy=0.0004, xcur=0.8, ylabel="Time (s)", xticks=[]);
tight_layout(); savefig("/Users/wenyue/Desktop/spt.pdf");


# function MultiStepForward_old(irz::Array{Int64,1}, irx::Array{Int64,1}, src::Source, fidMtx::FidMtx; tmax=0.5)
#     nz  =  src.nz;  nx = src.nx    ;
#     ext = src.ext;  iflag=src.iflag; dt = src.dt;
#     stl = src.ot ;  stu = src.ot+(src.nt-1)*dt;
#     nt  = round(Int64, tmax/dt)+1
#     spt1 = InitSnapShot(nz, nx, ext, iflag, dt, 1)
#     spt2 = InitSnapShot(nz, nx, ext, iflag, dt, 2)
#     AddSource!(spt1, src)
#     shotv = InitShotV(src.isz, src.isx, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
#     Spt2ShotV!(shotv, spt1)
#     for it = 2 : nt
#         OneStepForward!(spt2, spt1, fidMtx)
#         if stl <= (it-1)*dt <= stu
#            AddSource!(spt2, src)
#         end
#         CopySnapShot!(spt1, spt2)
#         Spt2ShotV!(shotv, spt1)
#     end
#     return shotv
# end
# @time shotv1 = MultiStepForward_old(irz, irx, src, fidMtx, tmax=tmax)
# @time shotv1 = MultiStepForward_old(irz, irx, src, fidMtx, tmax=tmax)

# SeisPlot(shotv.Vx, pclip=90)
#
# root = homedir(); path = join([root "/Desktop/wfd.bin"]);
# MultiStepForward(path, src, fidMtx, tmax=tmax);
# shotv1 = binWfd2shotv(irz, irx, path);
#
# path1 = join([root "/Desktop/moive"]);
# waveAnim(path1, path, cpt="Vx")
# fid = open(path, "r")
# lwfd = nz*nx*sizeof(Float32)*5; pre = sizeof(Float32)*5;
# nt = floor(Int64, (filesize(fid)-pre)/lwfd)
# shot = zeros(nt, length(irx))
# for it = 1 : nt
#     pos = pre + (it-1)*lwfd + 0*sizeof(Float32)*nz*nx
#     seek(fid, pos)
#     vx = reshape(convert(Array{Float64}, read(fid, Float32, nz*nx)), nz, nx)
#     SeisPlot(vx)
#     for ir = 1 : length(irx)
#         shot[it, ir] = vx[irz[ir], irx[ir]]
#     end
# end
#
# it = 300;
# wfd = ReadWfd(path, it)
# vx  = reshape(wfd.Vx, nz, nx);
# SeisPlot(vx)
# tmp = 0.0
# for ir = 1 : length(irx)
#     tmp += (shotv1.Vx[it, ir]-vx[irz[ir],irx[ir]])^2
# end

# path = "/Users/wenlei/Desktop/test.bin"
# MultiStepForward(path, src, fidMtx, dz, dx)

# it = 200;
# (vxx, vxz, vzx, vzz) = ReadPv(path, it);
#
# wfd = ReadWfd(path1, it)
# (vxx1, vxz1, vzx1, vzz1) = Dxz(wfd);

# wfd = ReadWfd(path, 300);
# d = wfd.Txx
# imshow(d, cmap="PuOr", vmax=maximum(abs(d)), vmin=-maximum(abs(d)));

# irx = collect(1:5:nx); irz = 5*ones(Int64, length(irx));
# shot = MultiStepForward(irz, irx, src, fidMtx);
