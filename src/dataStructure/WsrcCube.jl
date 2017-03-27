"""
nz, nx is the size of window, oz, ox is the starting index
temporarily is added as explosive source
"""
type WSrcCube
     oz :: Int64
     ox :: Int64
     nz :: Int64
     nx :: Int64
     ext:: Int64
     iflag :: Int64
     ot :: Float64
     dt :: Float64
     nt :: Int64
     p  :: Array{Float64, 3}
end

function InitWSrcCube(oz::Int64, ox::Int64, nz::Int64, nx::Int64, ext::Int64, iflag::Int64, ot::Float64, dt::Float64, nt::Int64)
    wsc = WSrcCube(oz, ox, nz, nx, ext, iflag, ot, dt, nt, zeros(nz, nx, nt))
    return wsc
end

function CopyWsc(wsc::WSrcCube)
   wsc1 = WSrcCube(wsc.oz, wsc.ox, wsc.nz, wsc.nx, wsc.ext, wsc.iflag, wsc.ot, wsc.dt, wsc.nt, copy(wsc.p))
   return wsc1
end

function l112norm(wsc::WSrcCube)
    (nz, nx, nt) = size(wsc.p)
    l = 0.0
    for ix = 1 : nx
        for iz = 1 : nz
            tmp = 0.0
            for it = 1 : nt
                tmp = tmp + wsc.p[iz,ix,it]*wsc.p[iz,ix,it]
            end
            l = l + sqrt(tmp)
        end
    end
    return l
end

function EngDisWsc(wsc::WSrcCube)
    (nz, nx, nt) = size(wsc.p)
    l = zeros(nz, nx)
    for ix = 1 : nx
        for iz = 1 : nz
            tmp = 0.0
            for it = 1 : nt
                tmp = tmp + wsc.p[iz,ix,it]*wsc.p[iz,ix,it]
            end
            l[iz,ix] = sqrt(tmp)
        end
    end
    return l
end

"""
convert explosive sources to WSrcCube
"""
function Srcs2Wsc(oz::Int64, ox::Int64, nz::Int64, nx::Int64, srcs::Array{Source,1})
    for is = 1 : length(srcs)
        if !srcs[is].Tzz_flag || !srcs[is].Txx_flag
           error("non explosive sources, can not convert srcs to WSrcCube")
        end
    end
    ext = srcs[1].ext; iflag = srcs[1].iflag;
    dt = srcs[1].dt
    (tl, tu) = SourcesTimeRange(srcs)
    nt = floor(Int64, (tu-tl)/dt) + 1
    wsc = InitWSrcCube(oz, ox, nz, nx, ext, iflag, tl, dt, nt)
    for it = 1 : nt
        for is = 1 : length(srcs)
            if srcs[is].ot <= tl+(it-1)*dt <=srcs[is].ot+(srcs[is].nt-1)*dt
               iz = srcs[is].isz; iz = iz - oz + 1;
               ix = srcs[is].isx; ix = ix - ox + 1;
               indt = it + round(Int64, (tl-srcs[is].ot)/dt)
               wsc.p[iz, ix, it] = wsc.p[iz, ix, it] + srcs[is].Tzz[indt]
            end
        end
    end
    return wsc
end

"""
sources is added as explosive sources
"""
function AddWsc2spt!(spt::SnapShot, wsc::WSrcCube)
    nz=spt.nz; nx=spt.nx; ext=spt.ext; iflag=spt.iflag;
    if iflag == 1
       zupper = ext
       Nz  = nz + 2*ext
    elseif iflag == 2
       zupper = 0
       Nz  = nz +   ext
    end
    Nx = nx + 2*ext
    it = spt.it
    dt = spt.dt
    indt = it - round(Int64, wsc.ot/dt)
    tmp  = zeros(Nz, Nx)
    izl  = zupper + wsc.oz; izu = zupper + wsc.oz + wsc.nz-1
    ixl  = ext    + wsc.ox; ixu = ext    + wsc.ox + wsc.nx-1
    tmp[izl:izu, ixl:ixu] = 1/2*wsc.p[:,:,indt]
    tmp = vec(tmp)
    spt.Tzzz[:] = spt.Tzzz[:] + tmp[:]
    spt.Tzzx[:] = spt.Tzzx[:] + tmp[:]
    spt.Txxz[:] = spt.Txxz[:] + tmp[:]
    spt.Txxx[:] = spt.Txxx[:] + tmp[:]
    return nothing
end

function spt2Wsc!(wsc::WSrcCube, spt::SnapShot)
    nz=spt.nz; nx=spt.nx; ext=spt.ext; iflag=spt.iflag;
    if iflag == 1
       zupper = ext
       Nz  = nz + 2*ext
    elseif iflag == 2
       zupper = 0
       Nz  = nz +   ext
    end
    Nx = nx + 2*ext
    it = spt.it
    dt = spt.dt
    indt = it - round(Int64, wsc.ot/dt)
    izl  = zupper + wsc.oz; izu = zupper + wsc.oz + wsc.nz-1
    ixl  = ext    + wsc.ox; ixu = ext    + wsc.ox + wsc.nx-1
    tmp  = 1/2*(spt.Tzzz + spt.Tzzx + spt.Txxz + spt.Txxx)
    tmp  = reshape(tmp, Nz, Nx)
    wsc.p[:, :, indt] = tmp[izl:izu, ixl:ixu]
    return nothing
end
