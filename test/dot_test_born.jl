using PyPlot, ElasticWave
nz = 177; nx = 333; ext= 10;   iflag = 1;
dx = 5. ; dz = 5. ; dt = 5e-4; tmax = 0.5; f0=30.;
vp =4000.*ones(nz,nx);  vs=2000.*ones(nz,nx); rho=2.5*ones(nz,nx);
fidMtx = CreateFidMtx(nz, nx, ext, iflag, vp, vs, rho, dz, dx, dt, f0);

root = homedir();
path_pv = join([root "/Desktop/pv.bin"])
isz = 1; isx = 165; ot = 0.0;
flags = vec([false false true true false]);
src = InitSource(isz, isx, nz, nx, ext, iflag, f0, ot, dt, flags);
MultiStepForward(path_pv, src, fidMtx, dz, dx)

path_wfd = join([root "/Desktop/wfd.bin"])
@time MultiStepForward(path_wfd, src, fidMtx, otype="wfd")

wfd = ReadWfd(path_wfd, 500); (vxx, vxz, vzx, vzz) = ReadPv(path_pv, 500);
tmp = reshape(wfd.Vx, nz, nx); SeisPlot(tmp)
tmp = reshape(wfd.Vz, nz, nx); SeisPlot(tmp)
SeisPlot(vzz+vxx); SeisPlot(vzx-vxz);

dm = zeros(nz, nx); du = zeros(nz, nx); dm[88, 110:220] = 500.; du[88, 110:220] = 500.;
irx = collect(1:5:nx); irz = 1*ones(Int64, length(irx));
@time shotv = MultiStepForward(irz, irx, path_pv,  dm, du, fidMtx; tmax = tmax)
path_wfd_born = join([root "/Desktop/wfd_born.bin"])
MultiStepForward(path_wfd_born, path_pv, dm, du, fidMtx; tmax=tmax);
wfd = ReadWfd(path_wfd_born, 300); SeisPlot(reshape(wfd.Vx, nz, nx)); SeisPlot(reshape(wfd.Vz, nz, nx));

SeisPlot(shotv.Vx); SeisPlot(shotv.Vz)

nt = shotv.nt
shotv_adj = InitShotV(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt);
wlet = Ricker(f0, dt); nw = length(wlet);
for ir = 1 : length(irz)
    if rand() < 0.9
       tmp = randn() * wlet
       inds = floor(Int64, rand()*(nt-nw))
       shotv_adj.Vx[inds:inds+nw-1 ,ir] = tmp[:]
    end
end
# imshow(shotv_adj.Vx, aspect=0.04);
(dm1, du1) = MultiStepAdjoint(shotv_adj, fidMtx, path_pv);
imshow(dm1); figure(); imshow(du1);

tmp1 = dot(vec(dm1), vec(dm)) + dot(vec(du1), vec(du))
tmp  = dot(vec(shotv.Vx), vec(shotv_adj.Vx)) + dot(vec(shotv.Vz), vec(shotv_adj.Vz))
# path_spt_adj = join([root "/Desktop/spt_adj.bin"])
# MultiStepAdjoint_test(path_spt_adj, shotv_adj, fidMtx)
# spt = ReadSnapShot(path_spt_adj, 500);
# imshow(reshape(spt.Vxx+spt.Vxz, 263, 337))
#
# dm1 = zeros(nz, nx); du1 = zeros(nz, nx);
# test1(dm1, du1, path_spt_adj, path_pv)
