using PyPlot, Seismic, ElasticWave

nz = 100; nx = 100; # model size
ext= 20; # thickness of PML layers
iflag = 1; #when iflag=1, four side obsorbing boundary, iflag=2, free surface
dx = 5. ; dz = 5. ; dt = 5e-4; f0=30.; tmax=0.5; # recording time and dominant frequency of source wavelet
rho=2.5e-4*ones(nz,nx);
# vp=2500.*ones(nz,nx);  vs=2500./sqrt(3)*ones(nz,nx);
# vp[20:end,:]=3000.; vs[20:end,:] = 3000./sqrt(3);
vp=3000.*ones(nz,nx);  vs=3000./sqrt(3)*ones(nz,nx);
vp[70:end,:]=4000.;  vs[70:end,:]=4000./sqrt(3); #physical model
fidMtx = CreateFidMtx(nz, nx, ext, iflag, vp, vs, rho, dz, dx, dt, f0); #finite difference matrix

# single source
isz = [50, 50]; isx = [35, 65]; ot = [0.0, 0.0]; #source location and starting time
flags = vec([false false true true false]); # add source to normal stress, simulate explosive sources
srcs = InitMultiSources(isz, isx, nz, nx, ext, iflag, f0, ot, dt, flags);

woz = 1; wox=1; wnz = 100; wnx = 100;
wsc = Srcs2Wsc(woz, wox, wnz, wnx, srcs);
wnt = size(wsc.p, 3)

irx = collect(1:2:nx); irz = ones(Int64,length(irx));
tmp = collect(1:2:nz); irz = vcat(irz, tmp); irx = vcat(irx, 5*ones(Int64,length(tmp)));
shot = MultiStepForward(nz, nx, irz, irx, wsc, fidMtx,  tmax=0.5);

# wsc1 = MultiStepAdjoint(woz, wox, wnz, wnx, wnt, shot, fidMtx);
shot.Vz = SeisAddNoise(shot.Vz, 0.1, L=5);
shot.Vx = SeisAddNoise(shot.Vx, 0.1, L=5);
SeisPlot(shot.Vz); SeisPlot(shot.Vx);
# wrand = RandWsc(woz, wox, wnz, wnx, ext, iflag, dt, wnt, f0);
# lambda = power_wsc(nz, nx, irz, irx, fidMtx, wrand, tmax; niter=30);

lambda = 2620.; mu=2.0;
(J, winv) = FISTA(woz, wox, wnz, wnx, wnt,
                  shot, fidMtx, mu, lambda; niter=30);
p = EngDisWsc(winv); SeisPlot(p, pclip=100);

shot1 = MultiStepForward(nz, nx, irz, irx, winv, fidMtx,  tmax=0.5);


function RandWsc(oz::Int64, ox::Int64, nz::Int64, nx::Int64,
                 ext::Int64, iflag::Int64, dt::Float64, nt::Int64,
                 f0::Float64)
    w = Ricker(f0, dt)
    l = floor(Int64, length(w)/2)
    a = randn(nz, nx, nt)
    for iz = 1 : nz
        for ix = 1 : nx
            tmp = vec(a[iz,ix,:])
            tmp1= conv(tmp, w)
            a[iz,ix,:] = tmp1[l+1:end-l]
        end
    end
    wsc = WSrcCube(oz, ox, nz, nx, ext, iflag, 0.0, dt, nt, a)
    return wsc
end

function power_wsc(nz::Int64, nx::Int64, irz::Array{Int64,1}, irx::Array{Int64,1},
                   fidMtx::FidMtx, wsc::WSrcCube, tmax::Float64; niter=20)
    # compute the dominant eigenvalue of A'A
    woz = wsc.oz; wox = wsc.ox; wnz = wsc.nz; wnx = wsc.nx; wnt =wsc.nt;
    lambda = 0.0
    for iter = 1 :niter
        shot = MultiStepForward(nz, nx, irz, irx, wsc, fidMtx, tmax=tmax)
        # shot.p = shot.p / sqrt(680)
        wsc  = MultiStepAdjoint(woz, wox, wnz, wnx, wnt, shot, fidMtx)
        # wsc.p = wsc.p / sqrt(680)
        lambda = vecnorm(wsc.p)
        wsc.p  = wsc.p / lambda
        println("iteration: $iter, maximum eig: $lambda")
    end
    return lambda
end

function softThresh_GP!(wsc::WSrcCube, T::Float64)
    (nz, nx, nt) = size(wsc.p)
    w = zeros(nz, nx)
    for ix = 1 : nx
        for iz = 1 : nz
            for it = 1 : nt
                w[iz, ix] = w[iz,ix] + wsc.p[iz,ix,it] * wsc.p[iz,ix,it]
            end
            w[iz,ix] = sqrt(w[iz,ix])
            if w[iz,ix] < T
               w[iz,ix] = 0.
            else
               w[iz,ix] = (w[iz,ix]-T) / w[iz,ix]
            end
        end
    end
    for it = 1 : nt
        wsc.p[:,:,it] = wsc.p[:,:,it] .* w
    end
    return nothing
end

function FISTA(woz::Int64, wox::Int64, wnz::Int64, wnx::Int64, wnt::Int64,
               d::ShotV, fidMtx::FidMtx, mu::Float64, lambda::Float64; niter=20)
    J = zeros(niter + 1)
    T = mu / (2*lambda)
    nz = d.nz; nx = d.nx;
    x0 = InitWSrcCube(woz, wox, wnz, wnx, d.ext, d.iflag, 0.0, d.dt, wnt)
    dtmp = InitShotV(0, 0, d.nz, d.nx, d.ext, d.iflag, d.irz, d.irx, 0.0, d.dt, d.nt)
    xk = CopyWsc(x0)
    yk = CopyWsc(x0)
    t0 = 1.0
    irz = d.irz; irx = d.irx; tmax=(d.nt-1)*d.dt
    for k = 1 : niter
        Gy = MultiStepForward(nz, nx, irz, irx, yk, fidMtx, tmax=tmax)
        dtmp.Vz = shot.Vz - Gy.Vz;   dtmp.Vx = shot.Vx - Gy.Vx;
        cost  = (vecnorm(dtmp.Vz))^2 + (vecnorm(dtmp.Vx))^2 + l112norm(yk);
        println("iterations: $k, cost $cost"); J[k] = cost;
        mtmp = MultiStepAdjoint(woz, wox, wnz, wnx, wnt, dtmp, fidMtx)
        xk.p = yk.p + mtmp.p/lambda
        softThresh_GP!(xk, T)
        t = (1 + sqrt(1+4*t0^2)) / 2
        yk.p = xk.p + (t0-1)/t*(xk.p-x0.p)
        t0 = t
        x0.p = copy(xk.p)
    end
    return J, xk
end

# dot product test
# forward
wsc.p = randn(size(wsc.p));
shot  = MultiStepForward(nz, nx, irz, irx, wsc, fidMtx,  tmax=0.5);

# adjoint
shot1 = InitShotV(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, size(shot.Vx,1));
shot1.Vx = randn(size(shot1.Vx));
shot1.Vz = randn(size(shot1.Vz));

wsc1  = MultiStepAdjoint(woz, wox, wnz, wnx, wnt, shot1, fidMtx);

# dot product
tmp = dot(vec(wsc.p), vec(wsc1.p))
tmp1= dot(vec(shot.Vx), vec(shot1.Vx)) + dot(vec(shot.Vz), vec(shot1.Vz))
