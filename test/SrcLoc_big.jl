using PyPlot, Seismic, ElasticWave

nz = 200; nx = 200; # model size
ext= 20; # thickness of PML layers
iflag = 1; #when iflag=1, four side obsorbing boundary, iflag=2, free surface
dx = 5. ; dz = 5. ; dt = 5e-4; f0=30.; tmax=0.8; # recording time and dominant frequency of source wavelet
rho=2.5e-4*ones(nz,nx);
# vp=2500.*ones(nz,nx);  vs=2500./sqrt(3)*ones(nz,nx);
# vp[20:end,:]=3000.; vs[20:end,:] = 3000./sqrt(3);
vp=3000.*ones(nz,nx);  vs=3000./sqrt(3)*ones(nz,nx);
vp[150:end,:]=4000.;  vs[150:end,:]=4000./sqrt(3); #physical model
fidMtx = CreateFidMtx(nz, nx, ext, iflag, vp, vs, rho, dz, dx, dt, f0); #finite difference matrix

# single source
isz = [100, 100]; isx = [85, 115]; ot = [0.0, 0.1]; #source location and starting time
flags = vec([false false true true false]); # add source to normal stress, simulate explosive sources
srcs = InitMultiSources(isz, isx, nz, nx, ext, iflag, f0, ot, dt, flags);

woz = 51; wox=51; wnz = 100; wnx = 100;
wsc = Srcs2Wsc(woz, wox, wnz, wnx, srcs);
wnt = 500; tmp = zeros(wnt,2);
tmp[1:346,1] = wsc.p[50,35,:];
tmp[1:346,2] = wsc.p[50,65,:];

fig = figure("1", figsize=(12,6)); s=25;
subplot(1,2,1); Seismic.SeisPlot(fignum="1", vp, cmap="jet", vmax=4000, vmin=3000, dx=5, dy=5, xlabel="X (m)", ylabel="Y (m)", ticksize=s, labelsize=s);
ax = gca(); ax[:set_xticks]([200, 400, 600, 800]); ax[:set_yticks]([200, 400, 600, 800]);
ax[:text](-100,60, "a)", fontsize=s, fontweight="bold")
subplot(1,2,2); Seismic.SeisPlot(fignum="1", tmp, style="wiggles", xcur=1.0, dy=0.0005, xlabel="Trace number", ylabel="Time (s)", ticksize=s, labelsize=s);
ax = gca(); ax[:set_xticks]([0, 1]); ax[:set_yticks]([0.05, 0.10, 0.15, 0.20]);
ax[:text](-1.3,0.015, "b)", fontsize=s, fontweight="bold")
tight_layout()

irx = collect(1:4:nx); irz = ones(Int64,length(irx));
tmp = collect(1:4:nz); irz = vcat(irz, tmp); irx = vcat(irx, 5*ones(Int64,length(tmp)));
shot = MultiStepForward(nz, nx, irz, irx, wsc, fidMtx,  tmax=tmax);
tmpx = shot.Vx; tmpz = shot.Vz;
shot.Vz = SeisAddNoise(shot.Vz, 0.5, L=5);
shot.Vx = SeisAddNoise(shot.Vx, 0.5, L=5);

fig = figure("2", figsize=(12,12)); s=15;
subplot(2,2,1); SeisPlot(fignum="2", tmpx, cmap="gray", pclip=95, dy=0.0005, ylabel="Time (s)", ticksize=s, labelsize=s)
ax = gca(); ax[:set_xticks]([20, 40, 60, 80]); ax[:set_yticks]([0.2, 0.4, 0.6]);
ax[:text](-8,0.04, "a)", fontsize=20, fontweight="bold")

subplot(2,2,2); SeisPlot(fignum="2", tmpz, cmap="gray", pclip=95, dy=0.0005, ticksize=s, labelsize=s)
ax = gca(); ax[:set_xticks]([20, 40, 60, 80]); ax[:set_yticks]([0.2, 0.4, 0.6]);
ax[:text](-8,0.04, "b)", fontsize=20, fontweight="bold")

subplot(2,2,3); SeisPlot(fignum="2", shot.Vx, cmap="gray", pclip=95, dy=0.0005, ylabel="Time (s)", xlabel="Trace number", ticksize=s, labelsize=s)
ax = gca(); ax[:set_xticks]([20, 40, 60, 80]); ax[:set_yticks]([0.2, 0.4, 0.6]);
ax[:text](-8,0.04, "c)", fontsize=20, fontweight="bold")

subplot(2,2,4); SeisPlot(fignum="2", shot.Vz, cmap="gray", pclip=95, dy=0.0005, xlabel="Trace number", ticksize=s, labelsize=s)
ax = gca(); ax[:set_xticks]([20, 40, 60, 80]); ax[:set_yticks]([0.2, 0.4, 0.6]);
ax[:text](-8,0.04, "d)", fontsize=20, fontweight="bold")

# wrand = RandWsc(woz, wox, wnz, wnx, ext, iflag, dt, wnt, f0);
# lambda = power_wsc(nz, nx, irz, irx, fidMtx, wrand, tmax; niter=30);

lambda = 1415.; mu=0.1;
(J, winv) = FISTA(woz, wox, wnz, wnx, wnt, shot, fidMtx, mu, lambda; niter=50);
lambda = 1415.; mu=0.01;
(J1, winv1) = FISTA(woz, wox, wnz, wnx, wnt, shot, fidMtx, mu, lambda; niter=50);
wadj = MultiStepAdjoint(1, 1, 200, 200, wnt, shot, fidMtx);
a = EngDisWsc(wadj);
p = EngDisWsc(winv); p1 = zeros(nz,nx); p1[woz:woz+wnz-1, wox:wox+wnx-1] = p[:,:];
v = maximum(p1)
SeisPlot(p1, vmax=v, vmin=0, cmap="jet", pclip=100)
v1 = maximum(a)
SeisPlot(a, vmax=v1, vmin=0, cmap="jet", pclip=100)
