using PyPlot, ElasticWave

# test add born source and its adjoint
path_pv  = "/Users/wenlei/Desktop/pv.bin"
path_spt = "/Users/wenlei/Desktop/spt.bin"
function ConvertBornSource2Spts(path_spt::String, path_pv::String, dm::Array{Float64,2}, du::Array{Float64,2})
    (nz, nx, ext, iflag, dt, ns) = InfoPv(path_pv)
    if iflag == 1
       Nz = nz + 2*ext
       upper = ext
     elseif iflag == 2
       Nz = nz +   ext
       upper = 0
    end
    Nx = nx + 2*ext
    # add source to Txx
    fid = open(path_spt, "w");
    write(fid, nz, nx, ext, iflag, dt)
    for it = 1 : ns
        (vxx, vxz, vzx, vzz) = ReadPv(path_pv, it)
        Txx = zeros(Nz, Nx)
        Txx[upper+1:nz+upper, ext+1:nx+ext] = 1/2*((vxx+vzz).*dm + 2*vxx.*du)
        Txx = vec(Txx)
        Tzz = zeros(Nz, Nx)
        Tzz[upper+1:nz+upper, ext+1:nx+ext] = 1/2*((vxx+vzz).*dm + 2*vzz.*du)
        Tzz = vec(Tzz)
        Txz = zeros(Nz, Nx)
        Txz[upper+1:nz+upper, ext+1:nx+ext] = 1/2*((vxz+vzx).*du)
        Txz = vec(Txz)
        write(fid, zeros(Nz*Nx))
        write(fid, zeros(Nz*Nx))
        write(fid, zeros(Nz*Nx))
        write(fid, zeros(Nz*Nx))
        write(fid, Txx         )
        write(fid, Txx         )
        write(fid, Tzz         )
        write(fid, Tzz         )
        write(fid, Txz         )
        write(fid, Txz         )
    end
    close(fid)
    return nothing
end

dm = 100*randn(nz, nx); du = 100*randn(nz, nx);
ConvertBornSource2Spts(path_spt, path_pv, dm, du);

function makeRandSnapShots(path_spt, nz, nx, ext, iflag, dt, nt)
    if iflag == 1
       Nz = nz + 2*ext
       upper = ext
     elseif iflag == 2
       Nz = nz +   ext
       upper = 0
    end
    Nx = nx + 2*ext
    fid = open(path_spt, "w")
    write(fid, nz, nx, ext, iflag, dt);
    for it = 1 : ns
        write(fid, randn(Nz*Nx*10))
    end
    close(fid)
    return nothing
end
(nz, nx, ext, iflag, dt, ns) = InfoPv(path_pv)
path_spt1 = "/Users/wenlei/Desktop/spt1.bin"
makeRandSnapShots(path_spt1, nz, nx, ext, iflag, dt, ns)


(dm1, du1) = spts2Born(path_spt1, path_pv)
tmp = dot(vec(dm), vec(dm1)) + dot(vec(du), vec(du1))
tmp1 = IptSpts(path_spt, path_spt1)

# =======================================================

nz = 233; nx = 277; ext= 30;   iflag = 2;
dx = 5. ; dz = 5. ; dt = 5e-4; tmax = 0.5; f0=30.;
vp =4000.*ones(nz,nx);  vs=2000.*ones(nz,nx); rho=2.5*ones(nz,nx);
fidMtx = CreateFidMtx(nz, nx, ext, iflag, vp, vs, rho, dz, dx, dt, f0);

path_pv  = "/Users/wenlei/Desktop/pv.bin"
path_spt = "/Users/wenlei/Desktop/spt.bin"
dm = 100*randn(nz, nx); du = 100*randn(nz, nx);
ConvertBornSource2Spts(path_spt, path_pv, dm, du);
irx = collect(1:5:nx); irz = 2*ones(Int64, length(irx));
shotv = MultiStepForward(irz, irx, path_spt, fidMtx)
imshow(shotv.Vx, aspect=0.04)

nt = shotv.nt
shotv_adj = InitShotV(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
wlet = Ricker(f0, dt); nw = length(wlet)
for ir = 1 : length(irz)
    if rand() < 0.6
       tmp = randn() * wlet
       inds = floor(Int64, rand()*(nt-nw))
       shotv_adj.Vx[inds:inds+nw-1 ,ir] = tmp[:]
    end
end
imshow(shotv_adj.Vx, aspect=0.04);

path_spt_adj = "/Users/wenlei/Desktop/spt_adj.bin"
MultiStepAdjoint(path_spt_adj, shotv_adj, fidMtx)
ReverseOrderSnapShots(path_spt_adj)

tmp1 = IptSpts(path_spt_adj, path_spt)
tmp  = dot(vec(shotv.Vx), vec(shotv_adj.Vx)) + dot(vec(shotv.Vz), vec(shotv_adj.Vz))
