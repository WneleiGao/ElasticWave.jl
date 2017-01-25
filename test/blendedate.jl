using PyPlot, ElasticWave
nz = 133; nx = 203; ext= 20; iflag = 1;
dx = 5. ; dz = 5. ; dt = 5e-4; tmax = 0.5; f0=30.;
vp =4000.*ones(nz,nx);  vs=2000.*ones(nz,nx); rho=2.5*ones(nz,nx)*1e-4;
fidMtx = CreateFidMtx(nz, nx, ext, iflag, vp, vs, rho, dz, dx, dt, f0);

if iflag == 1
   Nz = nz + 2*ext
elseif iflag == 2
   Nz = nz +   ext
end
Nx = nx + 2*ext

root = homedir();
path_pv = join([root "/Desktop/pv.bin"])
isx = collect(1:2:nx); isz = ones(Int64,length(isx)); ot = zeros(length(isx));
flags = vec([false false true true false]);
srcs  = InitMultiSources(isz, isx, nz, nx, ext, iflag, f0, ot, dt, flags);
MultiStepForward(path_pv, srcs, fidMtx, dz, dx)
(vxx, vxz, vzx, vzz) = ReadPv(path_pv, 300); p = vxx+vzz; s=vxz-vzx;
SeisPlot(reshape(p, Nz, Nx)); SeisPlot(reshape(s, Nz, Nx));

path_wfd = join([root "/Desktop/wfd.bin"])
MultiStepForward(path_wfd, srcs, fidMtx)
wfd = ReadWfd(path_wfd, 300);
SeisPlot(reshape(wfd.Vx, Nz, Nx))

dm = zeros(Nz, Nx); du = zeros(Nz, Nx); dm[floor(Int64, Nz/2), 1:end] = 500.; du[floor(Int64, Nz/2), 1:end] = 500.;
dm = vec(dm); du = vec(du);
irx = collect(1:2:nx); irz = 2*ones(Int64, length(irx));
shotv = MultiStepForward(irz, irx, path_pv,  dm, du, fidMtx)
SeisPlot(shotv.Vz); SeisPlot(shotv.Vx);

(dm1, du1) = MultiStepAdjoint(shotv, fidMtx, path_pv);
SeisPlot(reshape(dm1, Nz, Nx)); SeisPlot(reshape(du1, Nz, Nx)); 

# ==============================================================================
using PyPlot, ElasticWave

function main()
    nz = 177; nx = 203; ext= 30;   iflag = 1;
    dx = 5. ; dz = 5. ; dt = 4e-4; tmax = 0.5; f0=30.;
    vp =3500.*ones(nz,nx);  vs=1750.*ones(nz,nx); rho=2.5*ones(nz,nx);
    vp[77:end,:] = 5000.; vs[77:end,:] = 2500.
    fidMtx = CreateFidMtx(nz, nx, ext, iflag, vp, vs, rho, dz, dx, dt, f0);

    isx = collect(1:5:nx); isz = ones(Int64,length(isx)); ot = zeros(length(isx));
    flags = vec([false false true true false]);
    srcs  = InitMultiSources(isz, isx, nz, nx, ext, iflag, f0, ot, dt, flags);

    irx   = collect(1:1:nx); irz = 1*ones(Int64, length(irx));
    shotv = MultiStepForward(irz, irx, srcs, fidMtx)
    imshow(shotv.Vx, aspect=0.2, vmax=maximum(abs(shotv.Vx)), vmin=-maximum(abs(shotv.Vx)), cmap="PuOr"); figure();
    imshow(shotv.Vz, aspect=0.2, vmax=maximum(abs(shotv.Vz)), vmin=-maximum(abs(shotv.Vz)), cmap="PuOr")
    root = homedir()
    path_shotv = join([root "/Desktop/shotv.bin"])
    WriteShotV(path_shotv, shotv);

    shotv_FA = removeFirstArrival();
    shotv.Vx = shotv.Vx - shotv_FA.Vx;
    shotv.Vz = shotv.Vz - shotv_FA.Vz;
    figure(); imshow(shotv.Vx, aspect=0.2, vmax=maximum(abs(shotv.Vx)), vmin=-maximum(abs(shotv.Vx)), cmap="PuOr");
    figure(); imshow(shotv.Vz, aspect=0.2, vmax=maximum(abs(shotv.Vz)), vmin=-maximum(abs(shotv.Vz)), cmap="PuOr");

    path_shotv_DFA = join([root "/Desktop/shotv_DFA.bin"])
    WriteShotV(path_shotv_DFA, shotv);

    path_pv = join([root "/Desktop/pv.bin"])
    MultiStepForward(path_pv, srcs, fidMtx, dz, dx)
    (dm, du) = MultiStepAdjoint(shotv, fidMtx, path_pv);
    figure(); a = quantile(abs(dm[:]),0.98); imshow(dm, vmax=a, vmin=-a, cmap="PuOr");
    figure(); a = quantile(abs(du[:]),0.98); imshow(du, vmax=a, vmin=-a, cmap="PuOr");

    vps = modSmooth(vp, 20); vss = modSmooth(vs,20);
    fidMtx_SM = CreateFidMtx(nz, nx, ext, iflag, vps, vss, rho, dz, dx, dt, f0);

    path_pvsm = join([root "/Desktop/pv_sm.bin"])
    MultiStepForward(path_pvsm, srcs, fidMtx_SM, dz, dx)
    (dm2, du2) = MultiStepAdjoint(shotv, fidMtx_SM, path_pvsm);
    figure(); a = quantile(abs(dm1[:]),0.98); imshow(dm1, vmax=a, vmin=-a, cmap="PuOr");
    figure(); a = quantile(abs(du1[:]),0.98); imshow(du1, vmax=a, vmin=-a, cmap="PuOr");

    figure(); a = quantile(abs(dm2[:]),0.98); imshow(dm2, vmax=a, vmin=-a, cmap="PuOr");
    figure(); a = quantile(abs(du2[:]),0.98); imshow(du2, vmax=a, vmin=-a, cmap="PuOr");
end

function removeFirstArrival()
    nz = 177; nx = 203; ext= 30;   iflag = 1;
    dx = 10. ; dz = 10. ; dt = 1e-3; tmax = 1.0; f0=30.;
    vp =3500.*ones(nz,nx);  vs=1750.*ones(nz,nx); rho=2.5*ones(nz,nx);
    fidMtx = CreateFidMtx(nz, nx, ext, iflag, vp, vs, rho, dz, dx, dt, f0);

    isx = collect(1:5:nx); isz = ones(Int64,length(isx)); ot = zeros(length(isx));
    flags = vec([false false true true false]);
    srcs  = InitMultiSources(isz, isx, nz, nx, ext, iflag, f0, ot, dt, flags);

    irx   = collect(1:1:nx); irz = 1*ones(Int64, length(irx));
    shotv = MultiStepForward(irz, irx, srcs, fidMtx)
    return shotv
end

(vxx, vxz, vzx, vzz) = ReadPv(path_pv, 200);
p = (vzz+vxx); s = (vxz-vzx);
figure(); a = quantile(abs(p[:]),0.98); imshow(p, vmax=a, vmin=-a, cmap="PuOr");
figure(); a = quantile(abs(s[:]),0.98); imshow(s, vmax=a, vmin=-a, cmap="PuOr");
