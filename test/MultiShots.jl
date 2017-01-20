# inject plane wave source
using PyPlot, ElasticWave
nz = 150; nx = 250; ext= 20  ; iflag = 1 ;
dx = 20.; dz = 20.; dt = 2e-3; tmax  = 4.; f0=10.;
vp = 2500.*ones(nz, nx)   ; vs = 1450.*ones(nz, nx); rho=2.5*ones(nz,nx);
vp[26:100, :] = 3000.; vs[26:100, :] = 1752.;
for i = 1 : 25
    vp[i+50, 75-i:75+i] = 3500.;
    vs[i+50, 175-i:175+i] = 2020.;
end
vp[101:125, :] = 4000.; vs[101:125, :] = 2309.;
vp[126:nz, :] = 4500.; vs[126:nz, :] = 2600.;
SeisPlot(vp, wbox=6.7, hbox=4, vmax=maximum(vp), vmin=minimum(vp), cmap="afmhot", dx=dx, dy=dz, name="/Users/wenyue/Desktop/vp.pdf")
SeisPlot(vs, wbox=6.7, hbox=4, vmax=maximum(vs), vmin=minimum(vs), cmap="afmhot", dx=dx, dy=dz, name="/Users/wenyue/Desktop/vs.pdf")
fidMtx = CreateFidMtx(nz, nx, ext, iflag, vp, vs, rho, dz, dx, dt, f0);

vp = 2500.*ones(nz, nx); vs = 1450.*ones(nz, nx); rho=2.5*ones(nz,nx);
fidMtx_dp = CreateFidMtx(nz, nx, ext, iflag, vp, vs, rho, dz, dx, dt, f0);

# inject plane wave, output records with removing direct wave
isx = collect(1:3:nx); isz = ones(Int64, length(isx));
flags = [false, false, true, true, false];
thetal = -pi/6; thetau = pi/6; np = 21; vps = vp[1];
irx = collect(1:nx); irz = ones(Int64, length(irx));
root = "/Users/wenyue/Desktop/plane_wave/"
p = collect(linspace(thetal, thetau, np));
tdind = zeros(Int64, length(isx))
ot    = zeros(Float64, length(isx))
for ip = 1 : np
    for is = 1 : length(isx)
        lx = 0.0
        if p[ip] > 0
           lx = (isx[is]-1)*dx
        else p[ip] < 0
             lx = (isx[end]-isx[is])*dx
        end
        tdind[is] = floor(Int64, abs(lx)*sin(abs(p[ip]))/vps/dt) + 1
        ot[is] = (tdind[is]-1)*dt
    end
    srcs = InitMultiSources(isz, isx, nz, nx, ext, iflag, f0, ot, dt, flags);
    shotv = MultiStepForward(irz, irx, srcs, fidMtx, tmax=tmax)
    shot_tmp = MultiStepForward(irz, irx, srcs, fidMtx_dp, tmax=tmax)
    shotv.Vx[:] = shotv.Vx - shot_tmp.Vx
    shotv.Vz[:] = shotv.Vz - shot_tmp.Vz
    path = join([root "shotv" "$ip" ".bin"])
    WriteShotV(path, shotv)
end

# produce source side wave field(vxx, vxz, vzx, vzz)
vps = modSmooth(vp, 20); vss = modSmooth(vs,20);
root1 = "/Users/wenyue/Desktop/plane_wave/SWFD/"
fidMtx_SM = CreateFidMtx(nz, nx, ext, iflag, vps, vss, rho, dz, dx, dt, f0);
p = collect(linspace(thetal, thetau, np));
tdind = zeros(Int64, length(isx))
ot    = zeros(Float64, length(isx))
for ip = 1 : np
    for is = 1 : length(isx)
        lx = 0.0
        if p[ip] > 0
           lx = (isx[is]-1)*dx
        else p[ip] < 0
             lx = (isx[end]-isx[is])*dx
        end
        tdind[is] = floor(Int64, abs(lx)*sin(abs(p[ip]))/vps[1]/dt) + 1
        ot[is] = (tdind[is]-1)*dt
    end
    srcs = InitMultiSources(isz, isx, nz, nx, ext, iflag, f0, ot, dt, flags);
    path1 = join([root1 "pv" "$ip" ".bin"])
    MultiStepForward(path1, srcs, fidMtx_SM, dz, dx, tmax=tmax)
end

root = "/Users/wenyue/Desktop/plane_wave/Shotv/"
root1= "/Users/wenyue/Desktop/plane_wave/SWFD/"
root2= "/Users/wenyue/Desktop/plane_wave/Image/"

for ip = 1 : np
    path = join([root "shotv" "$ip" ".bin"])
    path1= join([root1 "pv" "$ip" ".bin"])
    path2= join([root2 "image" "$ip" ".bin"])
    shotv = ReadShotV(path)
    (dm, du) = MultiStepAdjoint(shotv, fidMtx_SM, path1)
    # dm = laplaceFilter(dm)
    # du = laplaceFilter(du)
    WriteImage(path2, dm, du)
end


Ip = zeros(nz, nx);
Is = zeros(nz, nx);
for ip = 1 : np
    path2= join([root2 "image" "$ip" ".bin"])
    (dm, du) = ReadImage(path2)
    # dm = gec(dm); du = gec(du)
    # (dvp, dvs) = lambdamu2vpvs(dm, du, vps, vss, rho)
    Ip = Ip + dm
    Is = Is + du
end


























tmax = 4.; dt = 2e-3; p = -pi/4; xs = collect(1:3:nx); dx=20.; v1 = 2500.;
path = "/Users/wenyue/Desktop/shotv_s_-45.bin";
root = "/Users/wenyue/Desktop/shotv/shotv";
tdind = planeSourceWfd(tmax, path, root, p, xs, dx, vp)
ot = (tdind-1)*dt; isx=collect(1:3:200); isz=ones(isx);
srcs = InitMultiSources(isz, isx, nz, nx, ext, iflag, f0, ot, dt, flags);
path_pv = "/Users/wenyue/Desktop//pv_-45.bin"
vps = modSmooth(vp, 20); vss = modSmooth(vs,20);
fidMtx_SM = CreateFidMtx(nz, nx, ext, iflag, vps, vss, rho, dz, dx, dt, f0);
MultiStepForward(path_pv, srcs, fidMtx_SM, dz, dx, tmax=tmax)
shotv_S = ReadShotV(path)
(dm, du) = MultiStepAdjoint(shotv_S, fidMtx_SM, path_pv);
figure(); a = quantile(abs(dm[:]),0.98); imshow(dm, vmax=a, vmin=-a, cmap="PuOr");
figure(); a = quantile(abs(du[:]),0.98); imshow(du, vmax=a, vmin=-a, cmap="PuOr");



function planeSrcShotV(tmax::Float64, path::ASCIIString, root::ASCIIString, p::Float64, xs::Array{Int64,1}, dx::Float64, vp::Float64)
    pathin = join([root "1" ".bin"])
    (nr, nz, nx, ext, iflag, irz, irx, dt, nst) = InfoShotV(pathin)
    nt = ceil(Int64, tmax/dt) + 1
    shotv_S = InitShotV(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
    tdind = zeros(Int64, length(xs))
    for ix = 1 : length(xs)
        lx = 0.0
        if p > 0
           lx = (xs[ix]-1)*dx
        else p < 0
           lx = (xs[end]-xs[ix])*dx
        end
        tdind[ix] = floor(Int64, abs(lx)*sin(abs(p))/vp/dt) + 1
    end
    for ix = 1 : length(xs)
        vx = zeros(nt,nr); vz = zeros(nt,nr);
        pathin = join([root "$ix" ".bin"])
        shotv  = ReadShotV(pathin)
        if tdind[ix]+nst-1 > nt
           tdend = nt
           nsend = nt - tdind[ix] + 1
        else
           tdend = tdind[ix] + nst - 1
           nsend = nst
        end
        vx[tdind[ix]:tdend,1:end] = shotv.Vx[1:nsend,1:end];
        vz[tdind[ix]:tdend,1:end] = shotv.Vz[1:nsend,1:end];
        shotv_S.Vx = shotv_S.Vx + vx
        shotv_S.Vz = shotv_S.Vz + vz
        WriteShotV(path, shotv_S)
    end
    return tdind
end






irx   = collect(1:1:nx); irz = 1*ones(Int64, length(irx));
flags = vec([false false true true false]);
xs = collect(1:3:200); root = "/Users/wenlei/Desktop/shotv"
for is = 1 : length(xs)
    isx = xs[is]; isz = 1; ot = 0.0;
    src = InitSource(isz, isx, nz, nx, ext, iflag, f0, ot, dt, flags);
    shotv = MultiStepForward(irz, irx, src, fidMtx, tmax=tmax);
    shotv_dp = MultiStepForward(irz, irx, src, fidMtx_dp, tmax=tmax);
    shotv.Vx[:] =shotv.Vx - shotv_dp.Vx
    shotv.Vz[:] =shotv.Vz - shotv_dp.Vz
    path = join([root "$is" ".bin"])
    WriteShotV(path, shotv)
    println("finished $is shot")
end
path = "/Users/wenyue/Desktop/shotv/shotv34.bin"
# shotv = ReadShotV(path);
# tmp=quantile(abs(vec(shotv.Vx)),0.95); imshow(shotv.Vx, aspect=0.1, vmax=tmp, vmin=-tmp, cmap="PuOr"); figure();
# tmp=quantile(abs(vec(shotv.Vz)),0.95); imshow(shotv.Vz, aspect=0.1, vmax=tmp, vmin=-tmp, cmap="PuOr")




function planeSrcShotV(tmax::Float64, path::ASCIIString, root::ASCIIString, p::Float64, xs::Array{Int64,1}, dx::Float64, vp::Float64)
    pathin = join([root "1" ".bin"])
    (nr, nz, nx, ext, iflag, irz, irx, dt, nst) = InfoShotV(pathin)
    nt = ceil(Int64, tmax/dt) + 1
    shotv_S = InitShotV(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
    tdind = zeros(Int64, length(xs))
    for ix = 1 : length(xs)
        lx = 0.0
        if p > 0
           lx = (xs[ix]-1)*dx
        else p < 0
           lx = (xs[end]-xs[ix])*dx
        end
        tdind[ix] = floor(Int64, abs(lx)*sin(abs(p))/vp/dt) + 1
    end
    for ix = 1 : length(xs)
        vx = zeros(nt,nr); vz = zeros(nt,nr);
        pathin = join([root "$ix" ".bin"])
        shotv  = ReadShotV(pathin)
        if tdind[ix]+nst-1 > nt
           tdend = nt
           nsend = nt - tdind[ix] + 1
        else
           tdend = tdind[ix] + nst - 1
           nsend = nst
        end
        vx[tdind[ix]:tdend,1:end] = shotv.Vx[1:nsend,1:end];
        vz[tdind[ix]:tdend,1:end] = shotv.Vz[1:nsend,1:end];
        shotv_S.Vx = shotv_S.Vx + vx
        shotv_S.Vz = shotv_S.Vz + vz
        WriteShotV(path, shotv_S)
    end
    return tdind
end

path_shotv = "/Users/wenyue/Desktop/shotv/shotv33.bin";
path_pv    = "/Users/wenyue/Desktop/pv33.bin";
ot = 0.0; isx=97; isz=1;
src = InitSource(isz, isx, nz, nx, ext, iflag, f0, ot, dt, flags);
MultiStepForward(path_pv, src, fidMtx_SM, dz, dx, tmax=tmax)
shotv = ReadShotV(path_shotv)
(dm, du) = MultiStepAdjoint(shotv, fidMtx_SM, path_pv);
figure(); a = quantile(abs(dm[:]),0.98); imshow(dm, vmax=a, vmin=-a, cmap="PuOr");
figure(); a = quantile(abs(du[:]),0.98); imshow(du, vmax=a, vmin=-a, cmap="PuOr");
