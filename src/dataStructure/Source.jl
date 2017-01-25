type Source
     isz  :: Int64
     isx  :: Int64
     nz   :: Int64
     nx   :: Int64
     ext  :: Int64
     iflag:: Int64
     ot   :: Float64
     dt   :: Float64
     nt   :: Int64
     Vx_flag  :: Bool
     Vz_flag  :: Bool
     Txx_flag :: Bool
     Tzz_flag :: Bool
     Txz_flag :: Bool
     Vx  :: Array{Float64,1}
     Vz  :: Array{Float64,1}
     Txx :: Array{Float64,1}
     Tzz :: Array{Float64,1}
     Txz :: Array{Float64,1}
end

function Ricker(f0::Float64, dt::Float64)
	  nw = 2.2/f0/dt
	  nw = 2*floor(Int,nw/2)+1
	  nc = floor(Int,nw/2)
	  w  = zeros(nw)
	  k  = collect(1:nw)
    k  = vec(k)
	  alpha = (nc-k+1)*f0*dt*pi
	  beta = alpha.^2
	  w = ((1.-beta.*2).*exp(-beta))
	  return w
end

function InitSource(isz::Int64, isx::Int64, nz::Int64, nx::Int64, ext::Int64, iflag::Int64, f0::Float64, ot::Float64, dt::Float64, flags::Array{Bool,1})
    if iflag == 1
       Nz = nz + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    if flags[1]
       w_Vx = Ricker(f0, dt)
       nt = length(w_Vx)
    else
       w_Vx = zeros(1)
    end
    if flags[2]
       w_Vz = Ricker(f0, dt)
       nt = length(w_Vz)
    else
       w_Vz = zeros(1)
    end
    if flags[3]
       w_Txx = Ricker(f0, dt)
       nt = length(w_Txx)
    else
       w_Txx = zeros(1)
    end
    if flags[4]
       w_Tzz = Ricker(f0, dt)
       nt = length(w_Tzz)
    else
       w_Tzz = zeros(1)
    end
    if flags[5]
       w_Txz = Ricker(f0, dt)
       nt = length(w_Txz)
    else
       w_Txz = zeros(1)
    end
    src = Source(isz, isx, nz, nx, ext, iflag, ot, dt, nt,
                 flags[1], flags[2], flags[3], flags[4], flags[5],
                 w_Vx, w_Vz, w_Txx, w_Tzz, w_Txz)
    return src
end

function InitMultiSources(isz::Array{Int64,1}, isx::Array{Int64,1}, nz::Int64, nx::Int64, ext::Int64, iflag::Int64, f0::Float64, ot::Array{Float64,1}, dt::Float64, flags::Array{Bool,1})
    ns = length(isz)
    srcs = Array(Source, ns)
    for is = 1 : ns
        srcs[is] = InitSource(isz[is], isx[is], nz, nx, ext, iflag, f0, ot[is], dt, flags)
    end
    return srcs
end

function SrcRange(src::Source)
    tl = src.ot
    tu = src.ot + (src.nt-1)*src.dt
    return tl, tu
end

function SourcesTimeRange(srcs::Array{Source,1})
    tmax = 0.0
    tmin = 0.0
    ns = length(srcs)
    dt = srcs[1].dt
    for is = 1 : ns
        tmp = srcs[is].ot + (srcs[is].nt-1)*dt
        tmp1= srcs[is].ot
        if tmp > tmax
           tmax = tmp
        end
        if tmp1 < tmin
           tmin = tmp1
        end
    end
    return tmin, tmax
end

function AddSource!(spt::SnapShot, src::Source)
    isz= src.isz
    isx= src.isx
    nz = src.nz
    nx = src.nx
    ext= src.ext
    iflag = src.iflag
    dt = src.dt
    if iflag == 1
       Nz = nz + 2*ext
       upper = ext
    elseif iflag ==2
       Nz = nz +   ext
       upper = 0
    end
    Nx = nx + 2*ext
    it  = spt.it
    if src.ot <= (it-1)*dt <= src.ot+(src.nt-1)*dt
       indt = it - round(Int64, src.ot/dt)
       indz = src.isz + upper
       indx = src.isx + ext
       pos = (indx-1) * Nz + indz
       if src.Vx_flag
          spt.Vxx[pos] = spt.Vxx[pos] + src.Vx[indt]/2
          spt.Vxz[pos] = spt.Vxz[pos] + src.Vx[indt]/2
       end
       if src.Vz_flag
          spt.Vzx[pos] = spt.Vzx[pos] + src.Vz[indt]/2
          spt.Vzz[pos] = spt.Vzz[pos] + src.Vz[indt]/2
       end
       if src.Txx_flag
          spt.Txxx[pos] = spt.Txxx[pos] + src.Txx[indt]/2
          spt.Txxz[pos] = spt.Txxz[pos] + src.Txx[indt]/2
       end
       if src.Tzz_flag
          spt.Tzzx[pos] = spt.Tzzx[pos] + src.Tzz[indt]/2
          spt.Tzzz[pos] = spt.Tzzz[pos] + src.Tzz[indt]/2
       end
       if src.Txz_flag
          spt.Txzx[pos] = spt.Txzx[pos] + src.Txz[indt]/2
          spt.Txzz[pos] = spt.Txzz[pos] + src.Txz[indt]/2
       end
    end
    return nothing
end

function AddMultiSources!(spt::SnapShot, srcs::Array{Source,1})
    nz    = srcs[1].nz
    nx    = srcs[1].nx
    ext   = srcs[1].ext
    iflag = srcs[1].iflag
    dt = srcs[1].dt
    if iflag == 1
       Nz = nz + 2*ext
       upper = ext
    elseif iflag ==2
       Nz = nz +   ext
       upper = 0
    end
    Nx = nx + 2*ext
    it = spt.it
    ns = length(srcs)
    for is = 1: ns
        if srcs[is].ot <= (it-1)*dt <= srcs[is].ot+(srcs[is].nt-1)*dt
           indt = it - round(Int64, srcs[is].ot/dt)
           indz = srcs[is].isz + upper
           indx = srcs[is].isx + ext
           pos  = (indx-1) * Nz + indz
           if srcs[is].Vx_flag
              spt.Vxx[pos] = spt.Vxx[pos] + srcs[is].Vx[indt]/2
              spt.Vxz[pos] = spt.Vxz[pos] + srcs[is].Vx[indt]/2
           end
           if srcs[is].Vz_flag
              spt.Vzx[pos] = spt.Vzx[pos] + srcs[is].Vz[indt]/2
              spt.Vzz[pos] = spt.Vzz[pos] + srcs[is].Vz[indt]/2
           end
           if srcs[is].Txx_flag
              spt.Txxx[pos] = spt.Txxx[pos] + srcs[is].Txx[indt]/2
              spt.Txxz[pos] = spt.Txxz[pos] + srcs[is].Txx[indt]/2
           end
           if srcs[is].Tzz_flag
              spt.Tzzx[pos] = spt.Tzzx[pos] + srcs[is].Tzz[indt]/2
              spt.Tzzz[pos] = spt.Tzzz[pos] + srcs[is].Tzz[indt]/2
           end
           if srcs[is].Txz_flag
              spt.Txzx[pos] = spt.Txzx[pos] + srcs[is].Txz[indt]/2
              spt.Txzz[pos] = spt.Txzz[pos] + srcs[is].Txz[indt]/2
           end
        end
    end
    return nothing
end

function srcs2spt(path::String, srcs::Array{Source,1})
    (tl, tu) = SourcesTimeRange(srcs)
    nz = srcs[1].nz; nx = srcs[1].nx; ext = srcs[1].ext; iflag = srcs[1].iflag;
    dt = srcs[1].dt; nt = floor(Int64, tu/dt)+1;
    spt = InitSnapShot(nz, nx, ext, iflag, dt, 1)
    AddMultiSources!(spt, srcs)
    fid = WriteSnapShot(path, spt)
    for it = 2 : nt
        spt = InitSnapShot(nz, nx, ext, iflag, dt, it)
        AddMultiSources!(spt, srcs)
        WriteSnapShot(fid, spt)
    end
    close(fid)
    return nothing
end

function srcs2wfd(path::String, srcs::Array{Source,1})
    (tl, tu) = SourcesTimeRange(srcs)
    nz = srcs[1].nz; nx = srcs[1].nx; ext = srcs[1].ext; iflag = srcs[1].iflag;
    dt = srcs[1].dt; nt = floor(Int64, tu/dt)+1;
    wfd = InitWfd(nz, nx, ext, iflag, dt, 1)
    AddSrcs2Wfd!(wfd, srcs)
    fid = WriteSnapShot(path, spt)
    for it = 2 : nt
        spt = InitSnapShot(nz, nx, ext, iflag, dt, it)
        AddMultiSources!(spt, srcs)
        WriteSnapShot(fid, spt)
    end
    close(fid)
    return nothing
end

function AddSourceBorn!(spt::SnapShot, path::String, dm::Array{Float64,1}, du::Array{Float64,1})
    it = spt.it;
    nz = spt.nz; nx = spt.nx; ext = spt.ext;
    if spt.iflag == 1
       Nz = nz + 2*ext
       upper = ext
     elseif spt.iflag == 2
       Nz = nz +   ext
       upper = 0
    end
    Nx = nx + 2*ext
    (vxx, vxz, vzx, vzz) = ReadPv(path, it)
    pc = vxx+vzz
    sc = vxz+vzx
    # add source to Txx
    tmp = 1/2*(pc.*dm+2*vxx.*du);
    spt.Txxx = spt.Txxx + tmp
    spt.Txxz = spt.Txxz + tmp
    # add source to Tzz
    tmp = 1/2*(pc.*dm+2*vzz.*du)
    spt.Tzzx = spt.Tzzx + tmp
    spt.Tzzz = spt.Tzzz + tmp
    # add source to Txz
    tmp = 1/2*(sc.*du)
    spt.Txzx = spt.Txzx + tmp
    spt.Txzz = spt.Txzz + tmp
    return nothing
end
# function AddSourceBorn!(spt::SnapShot, path::String, dm::Array{Float64,2}, du::Array{Float64,2})
#     it = spt.it;
#     nz = spt.nz; nx = spt.nx; ext = spt.ext;
#     if spt.iflag == 1
#        Nz = nz + 2*ext
#        upper = ext
#      elseif spt.iflag == 2
#        Nz = nz +   ext
#        upper = 0
#     end
#     Nx = nx + 2*ext
#     # add source to Txx
#     tmp = zeros(Nz, Nx)
#     (vxx, vxz, vzx, vzz) = ReadPv(path, it)
#     tmp[upper+1:nz+upper, ext+1:nx+ext] = 1/2*((vxx+vzz).*dm+2*vxx.*du)
#     tmp = vec(tmp)
#     spt.Txxx = spt.Txxx + tmp
#     spt.Txxz = spt.Txxz + tmp
#     # add source to Tzz
#     tmp = zeros(Nz, Nx)
#     tmp[upper+1:nz+upper, ext+1:nx+ext] = 1/2*((vxx+vzz).*dm+2*vzz.*du)
#     tmp = vec(tmp)
#     spt.Tzzx = spt.Tzzx + tmp
#     spt.Tzzz = spt.Tzzz + tmp
#     # add source to Txz
#     tmp = zeros(Nz, Nx)
#     tmp[upper+1:nz+upper, ext+1:nx+ext] = 1/2*((vxz+vzx).*du)
#     tmp = vec(tmp)
#     spt.Txzx = spt.Txzx + tmp
#     spt.Txzz = spt.Txzz + tmp
#     return nothing
# end
function ConvertBornSource2Spts(path_spt::String, path_pv::String, dm::Array{Float64,1}, du::Array{Float64,1})
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
        pc = vxx+vzz
        Txx = 1/2*(pc.*dm + 2*vxx.*du)
        Tzz = 1/2*(pc.*dm + 2*vzz.*du)
        Txz = 1/2*(sc.*du)
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
# function ConvertBornSource2Spts(path_spt::String, path_pv::String, dm::Array{Float64,2}, du::Array{Float64,2})
#     (nz, nx, ext, iflag, dt, ns) = InfoPv(path_pv)
#     if iflag == 1
#        Nz = nz + 2*ext
#        upper = ext
#      elseif iflag == 2
#        Nz = nz +   ext
#        upper = 0
#     end
#     Nx = nx + 2*ext
#     # add source to Txx
#     fid = open(path_spt, "w");
#     write(fid, nz, nx, ext, iflag, dt)
#     for it = 1 : ns
#         (vxx, vxz, vzx, vzz) = ReadPv(path_pv, it)
#         Txx = zeros(Nz, Nx)
#         Txx[upper+1:nz+upper, ext+1:nx+ext] = 1/2*((vxx+vzz).*dm + 2*vxx.*du)
#         Txx = vec(Txx)
#         Tzz = zeros(Nz, Nx)
#         Tzz[upper+1:nz+upper, ext+1:nx+ext] = 1/2*((vxx+vzz).*dm + 2*vzz.*du)
#         Tzz = vec(Tzz)
#         Txz = zeros(Nz, Nx)
#         Txz[upper+1:nz+upper, ext+1:nx+ext] = 1/2*((vxz+vzx).*du)
#         Txz = vec(Txz)
#         write(fid, zeros(Nz*Nx))
#         write(fid, zeros(Nz*Nx))
#         write(fid, zeros(Nz*Nx))
#         write(fid, zeros(Nz*Nx))
#         write(fid, Txx         )
#         write(fid, Txx         )
#         write(fid, Tzz         )
#         write(fid, Tzz         )
#         write(fid, Txz         )
#         write(fid, Txz         )
#     end
#     close(fid)
#     return nothing
# end
