type SnapShot
     nz   :: Int64
     nx   :: Int64
     ext  :: Int64
     iflag:: Int64
     dt   :: Float64
     it   :: Int64
     Vxx  :: Array{Float64,1}
     Vxz  :: Array{Float64,1}
     Vzx  :: Array{Float64,1}
     Vzz  :: Array{Float64,1}
     Txxx :: Array{Float64,1}
     Txxz :: Array{Float64,1}
     Tzzx :: Array{Float64,1}
     Tzzz :: Array{Float64,1}
     Txzx :: Array{Float64,1}
     Txzz :: Array{Float64,1}
end

function InitSnapShot(nz::Int64, nx::Int64, ext::Int64, iflag::Int64, dt::Float64, it::Int64)
    if iflag == 1
       Nz = nz + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    spt= SnapShot(nz, nx, ext, iflag, dt, it,
                  zeros(Nz*Nx), zeros(Nz*Nx), zeros(Nz*Nx), zeros(Nz*Nx), zeros(Nz*Nx),
                  zeros(Nz*Nx), zeros(Nz*Nx), zeros(Nz*Nx), zeros(Nz*Nx), zeros(Nz*Nx))
    return spt
end

type Wfd
     nz   :: Int64
     nx   :: Int64
     ext  :: Int64
     iflag:: Int64
     dt   :: Float64
     it   :: Int64
     Vx   :: Array{Float64,1}
     Vz   :: Array{Float64,1}
     Txx  :: Array{Float64,1}
     Tzz  :: Array{Float64,1}
     Txz  :: Array{Float64,1}
end

function InitWfd(nz::Int64, nx::Int64, ext::Int64, iflag::Int64, dt::Float64, it::Int64)
    wfd= Wfd(nz, nx, ext, iflag, dt, it,
             zeros(nz*nx), zeros(nz*nx), zeros(nz*nx), zeros(nz*nx), zeros(nz*nx))
    return wfd
end

function CopySnapShot!(snapShot1, snapShot2)
    if snapShot2.nz != snapShot1.nz || snapShot2.nx != snapShot1.nx || snapShot2.ext != snapShot1.ext || snapShot2.iflag != snapShot1.iflag
       error("the two snapShot are different")
    end
    snapShot1.it   = snapShot2.it
    n = length(snapShot1.Vxx)
    @inbounds for i = 1 : n
        snapShot1.Vxx[i]  = snapShot2.Vxx[i]
    end
    @inbounds for i = 1 : n
        snapShot1.Vxz[i]  = snapShot2.Vxz[i]
    end
    @inbounds for i = 1 : n
        snapShot1.Vzx[i]  = snapShot2.Vzx[i]
    end
    @inbounds for i = 1 : n
        snapShot1.Vzz[i]  = snapShot2.Vzz[i]
    end
    @inbounds for i = 1 : n
        snapShot1.Txxx[i]  = snapShot2.Txxx[i]
    end
    @inbounds for i = 1 : n
        snapShot1.Txxz[i]  = snapShot2.Txxz[i]
    end
    @inbounds for i = 1 : n
        snapShot1.Tzzx[i]  = snapShot2.Tzzx[i]
    end
    @inbounds for i = 1 : n
        snapShot1.Tzzz[i]  = snapShot2.Tzzz[i]
    end
    @inbounds for i = 1 : n
        snapShot1.Txzx[i]  = snapShot2.Txzx[i]
    end
    @inbounds for i = 1 : n
        snapShot1.Txzz[i]  = snapShot2.Txzz[i]
    end
    return nothing
end

function CopySnapShot(snapShot::SnapShot)
    spt = SnapShot(snapShot.nz, snapShot.nx, snapShot.ext, snapShot.iflag, snapShot.dt, snapShot.it,
                   snapShot.Vxx , snapShot.Vxz , snapShot.Vzx , snapShot.Vzz,
                   snapShot.Txxx, snapShot.Txxz, snapShot.Tzzx, snapShot.Tzzz, snapShot.Txzx, snapShot.Txzz)
    return spt
end

function Add2SnapShots!(snapShot2::SnapShot, snapShot1::SnapShot)
    snapShot2.Vxx = snapShot2.Vxx + snapShot1.Vxx
    snapShot2.Vxz = snapShot2.Vxz + snapShot1.Vxz
    snapShot2.Vzx = snapShot2.Vzx + snapShot1.Vzx
    snapShot2.Vzz = snapShot2.Vzz + snapShot1.Vzz
    snapShot2.Txxx = snapShot2.Txxx + snapShot1.Txxx
    snapShot2.Txxz = snapShot2.Txxz + snapShot1.Txxz
    snapShot2.Tzzx = snapShot2.Tzzx + snapShot1.Tzzx
    snapShot2.Tzzz = snapShot2.Tzzz + snapShot1.Tzzz
    snapShot2.Txzx = snapShot2.Txzx + snapShot1.Txzx
    snapShot2.Txzz = snapShot2.Txzz + snapShot1.Txzz
    return nothing
end

function WriteSnapShot(pathout::String, snapShot::SnapShot)
    fid = open(pathout, "w")
    write(fid, convert(Int32  , snapShot.nz))
    write(fid, convert(Int32  , snapShot.nx))
    write(fid, convert(Int32  , snapShot.ext))
    write(fid, convert(Int32  , snapShot.iflag))
    write(fid, convert(Float32, snapShot.dt))
    write(fid, convert(Array{Float32}, snapShot.Vxx))
    write(fid, convert(Array{Float32}, snapShot.Vxz))
    write(fid, convert(Array{Float32}, snapShot.Vzx))
    write(fid, convert(Array{Float32}, snapShot.Vzz))
    write(fid, convert(Array{Float32}, snapShot.Txxx))
    write(fid, convert(Array{Float32}, snapShot.Txxz))
    write(fid, convert(Array{Float32}, snapShot.Tzzx))
    write(fid, convert(Array{Float32}, snapShot.Tzzz))
    write(fid, convert(Array{Float32}, snapShot.Txzx))
    write(fid, convert(Array{Float32}, snapShot.Txzz))
    return fid
end

function WriteSnapShot(fid::IOStream, snapShot::SnapShot)
    write(fid, convert(Array{Float32}, snapShot.Vxx))
    write(fid, convert(Array{Float32}, snapShot.Vxz))
    write(fid, convert(Array{Float32}, snapShot.Vzx))
    write(fid, convert(Array{Float32}, snapShot.Vzz))
    write(fid, convert(Array{Float32}, snapShot.Txxx))
    write(fid, convert(Array{Float32}, snapShot.Txxz))
    write(fid, convert(Array{Float32}, snapShot.Tzzx))
    write(fid, convert(Array{Float32}, snapShot.Tzzz))
    write(fid, convert(Array{Float32}, snapShot.Txzx))
    write(fid, convert(Array{Float32}, snapShot.Txzz))
    flush(fid)
    return nothing
end

function ReverseOrderSnapShots(pathin::String)
    fid = open(pathin, "r")
    pathin_tmp = join([pathin "_tmp"])
    fid1  = open(pathin_tmp, "w")
    nz    = read(fid, Int32)
    nx    = read(fid, Int32)
    ext   = read(fid, Int32)
    iflag = read(fid, Int32)
    dt    = read(fid, Float32)
    write(fid1, nz, nx, ext, iflag, dt);
    if iflag == 1
       Nz = nz + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    snapSize = sizeof(Float32)*Nz*Nx*10
    nt = round(Int64, ((filesize(fid)-sizeof(Float32)*5)/snapSize))
    for it = nt : -1 : 1
        position = sizeof(Float32)*5 + (it-1)*snapSize
        seek(fid, position)
        tmp = read(fid, Float32, Nz*Nx*10)
        write(fid1, tmp)
    end
    close(fid); close(fid1)
    mv(pathin_tmp, pathin, remove_destination=true)
    return nothing
end

function WriteWfd(pathout::String, snapShot::SnapShot)
    fid = open(pathout, "w")
    nz = snapShot.nz;  nx = snapShot.nx;
    ext= snapShot.ext; iflag = snapShot.iflag;
    write(fid, convert(Int32,  nz))
    write(fid, convert(Int32,  nx))
    write(fid, convert(Int32,  ext))
    write(fid, convert(Int32,  iflag))
    write(fid, convert(Float32, snapShot.dt))
    if iflag == 1
       Nz = nz + 2*ext
       zupper = ext
    elseif iflag == 2
       Nz = nz +   ext
       zupper = 0
    end
    Nx = nx + 2*ext
    tmp = snapShot.Vxx+snapShot.Vxz;   tmp = reshape(tmp, Nz, Nx); write(fid, convert(Array{Float32}, vec(tmp[zupper+1:zupper+nz, ext+1:ext+nx])))
    tmp = snapShot.Vzx+snapShot.Vzz;   tmp = reshape(tmp, Nz, Nx); write(fid, convert(Array{Float32}, vec(tmp[zupper+1:zupper+nz, ext+1:ext+nx])))
    tmp = snapShot.Txxx+snapShot.Txxz; tmp = reshape(tmp, Nz, Nx); write(fid, convert(Array{Float32}, vec(tmp[zupper+1:zupper+nz, ext+1:ext+nx])))
    tmp = snapShot.Tzzx+snapShot.Tzzz; tmp = reshape(tmp, Nz, Nx); write(fid, convert(Array{Float32}, vec(tmp[zupper+1:zupper+nz, ext+1:ext+nx])))
    tmp = snapShot.Txzx+snapShot.Txzz; tmp = reshape(tmp, Nz, Nx); write(fid, convert(Array{Float32}, vec(tmp[zupper+1:zupper+nz, ext+1:ext+nx])))
    return fid
end

function WriteWfd(fid::IOStream, snapShot::SnapShot)
    nz = snapShot.nz;  nx = snapShot.nx;
    ext= snapShot.ext; iflag = snapShot.iflag;
    if iflag == 1
       Nz = nz + 2*ext
       zupper = ext
    elseif iflag == 2
       Nz = nz +   ext
       zupper = 0
    end
    Nx = nx + 2*ext
    tmp = snapShot.Vxx+snapShot.Vxz;   tmp = reshape(tmp, Nz, Nx); write(fid, convert(Array{Float32}, vec(tmp[zupper+1:zupper+nz, ext+1:ext+nx])))
    tmp = snapShot.Vzx+snapShot.Vzz;   tmp = reshape(tmp, Nz, Nx); write(fid, convert(Array{Float32}, vec(tmp[zupper+1:zupper+nz, ext+1:ext+nx])))
    tmp = snapShot.Txxx+snapShot.Txxz; tmp = reshape(tmp, Nz, Nx); write(fid, convert(Array{Float32}, vec(tmp[zupper+1:zupper+nz, ext+1:ext+nx])))
    tmp = snapShot.Tzzx+snapShot.Tzzz; tmp = reshape(tmp, Nz, Nx); write(fid, convert(Array{Float32}, vec(tmp[zupper+1:zupper+nz, ext+1:ext+nx])))
    tmp = snapShot.Txzx+snapShot.Txzz; tmp = reshape(tmp, Nz, Nx); write(fid, convert(Array{Float32}, vec(tmp[zupper+1:zupper+nz, ext+1:ext+nx])))
    flush(fid)
    return nothing
end

function ReadWfd(pathin::String, it::Int64)
    fid = open(pathin, "r")
    nz  = convert(Int64, read(fid, Int32))
    nx  = convert(Int64, read(fid, Int32))
    ext = convert(Int64, read(fid, Int32))
    iflag = convert(Int64, read(fid, Int32))
    dt  = convert(Float64, read(fid, Float32))
    wfd = InitWfd(nz, nx, ext, iflag, dt, it)
    wfdSize = sizeof(Float32)*nz*nx*5
    position = sizeof(Float32)*5 + (it-1)*wfdSize
    seek(fid, position)
    wfd.Vx  = convert(Array{Float64}, read(fid,Float32, nz*nx))
    wfd.Vz  = convert(Array{Float64}, read(fid,Float32, nz*nx))
    wfd.Txx = convert(Array{Float64}, read(fid,Float32, nz*nx))
    wfd.Tzz = convert(Array{Float64}, read(fid,Float32, nz*nx))
    wfd.Txz = convert(Array{Float64}, read(fid,Float32, nz*nx))
    close(fid)
    return wfd
end

function InfoWfd(pathin::String; flag=false)
    fid = open(pathin, "r")
    nz  = convert(Int64, read(fid, Int32))
    nx  = convert(Int64, read(fid, Int32))
    ext = convert(Int64, read(fid, Int32))
    iflag = convert(Int64, read(fid, Int32))
    dt  = convert(Float64, read(fid, Float32))
    wfdSize = sizeof(Float32)*nz*nx*5
    nt = round(Int64, (filesize(fid)-sizeof(Int32)*5)/wfdSize)
    if flag == true
       println("nz: $nz")
       println("nx: $nx")
       println("ext: $ext")
       println("iflag: $iflag")
       println("dt: $dt")
       println("nt: $nt")
       tmax = (nt-1)*dt
       println("tmax: $tmax")
       fSize = filesize(fid)/1024/1024/1024
       println("file size: $fSize G")
    end
    close(fid)
    return nz, nx, ext, iflag, dt, nt
end

function InfoSnapShots(pathin::String; flag=false)
    fid = open(pathin, "r")
    nz  = convert(Int64, read(fid, Int32))
    nx  = convert(Int64, read(fid, Int32))
    ext = convert(Int64, read(fid, Int32))
    iflag = convert(Int64, read(fid, Int32))
    dt  = convert(Float64, read(fid, Float32))
    if iflag == 1
       Nz = nz + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    snapSize = sizeof(Float32)*Nz*Nx*10
    nt = round(Int64, (filesize(fid)-sizeof(Int32)*5)/snapSize)
    if flag == true
       println("nz: $nz")
       println("nx: $nx")
       println("ext: $ext")
       println("iflag: $iflag")
       println("dt: $dt")
       println("nt: $nt")
       tmax = (nt-1)*dt
       println("tmax: $tmax")
       fSize = filesize(fid)/1024/1024/1024
       println("file size: $fSize G")
    end
    close(fid)
    return nz, nx, ext, iflag, dt, nt
end

function ReadSnapShot(pathin::String, it::Int64)
    fid = open(pathin, "r")
    nz  = convert(Int64, read(fid, Int32))
    nx  = convert(Int64, read(fid, Int32))
    ext = convert(Int64, read(fid, Int32))
    iflag = convert(Int64, read(fid, Int32))
    dt  = convert(Float64, read(fid, Float32))
    if iflag == 1
       Nz = nz + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    snapShot = InitSnapShot(nz, nx, ext, iflag, dt, it)
    snapSize = sizeof(Float32)*Nz*Nx*10
    position = sizeof(Float32)*5 + (it-1)*snapSize
    seek(fid, position)
    snapShot.Vxx  = convert(Array{Float64}, read(fid, Float32, Nz*Nx))
    snapShot.Vxz  = convert(Array{Float64}, read(fid, Float32, Nz*Nx))
    snapShot.Vzx  = convert(Array{Float64}, read(fid, Float32, Nz*Nx))
    snapShot.Vzz  = convert(Array{Float64}, read(fid, Float32, Nz*Nx))
    snapShot.Txxx = convert(Array{Float64}, read(fid, Float32, Nz*Nx))
    snapShot.Txxz = convert(Array{Float64}, read(fid, Float32, Nz*Nx))
    snapShot.Tzzx = convert(Array{Float64}, read(fid, Float32, Nz*Nx))
    snapShot.Tzzz = convert(Array{Float64}, read(fid, Float32, Nz*Nx))
    snapShot.Txzx = convert(Array{Float64}, read(fid, Float32, Nz*Nx))
    snapShot.Txzz = convert(Array{Float64}, read(fid, Float32, Nz*Nx))
    close(fid)
    return snapShot
end

function ExtractCpt(path::String, pathin::String; itype="wfd", cpt="Vx")
    fid = open(pathin, "r")
    nz    = read(fid, Int32)
    nx    = read(fid, Int32)
    ext   = read(fid, Int32)
    iflag = read(fid, Int32)
    dt    = read(fid, Float32)
    fido = open(path, "w")
    write(fido, nz, nx, ext, iflag, dt)
    if itype == "spt"
       if iflag == 1
          Nz = nz + 2*ext
          upper = ext
       elseif iflag == 2
          Nz = nz +   ext
          upper = 0
       end
       Nx = nx + 2*ext
       sptSize = sizeof(Float32)*Nz*Nx*10
       nt = round(Int64, ((filesize(fid)-sizeof(Float32)*5)/sptSize))
       inds = 0
       if cpt == "Vx"
          inds = sizeof(Float32)*5
       elseif cpt == "Vz"
          inds = sizeof(Float32)*5 + 2*Nz*Nx*sizeof(Float32)
       elseif cpt == "Txx"
          inds = sizeof(Float32)*5 + 4*Nz*Nx*sizeof(Float32)
       elseif cpt == "Txx"
          inds = sizeof(Float32)*5 + 6*Nz*Nx*sizeof(Float32)
       elseif cpt == "Txx"
          inds = sizeof(Float32)*5 + 8*Nz*Nx*sizeof(Float32)
       end
       for it = 1: nt
           position = (it-1)*sptSize + inds
           seek(fid, position)
           tmp = read(fid,Float32,Nz*Nx) + read(fid,Float32,Nz*Nx)
           tmp = reshape(tmp, Int64(Nz), Int64(Nx))
           write(fido, vec(tmp[upper+1:upper+nz, ext+1:ext+nx]))
       end
    elseif itype == "wfd"
       sptSize = sizeof(Float32)*nz*nx*5
       nt = round(Int64, ((filesize(fid)-sizeof(Float32)*5)/sptSize))
       inds = 0
       if cpt == "Vx"
         inds = sizeof(Float32)*5
       elseif cpt == "Vz"
         inds = sizeof(Float32)*5 +   nz*nx*sizeof(Float32)
       elseif cpt == "Txx"
         inds = sizeof(Float32)*5 + 2*nz*nx*sizeof(Float32)
       elseif cpt == "Txx"
         inds = sizeof(Float32)*5 + 2*nz*nx*sizeof(Float32)
       elseif cpt == "Txx"
         inds = sizeof(Float32)*5 + 2*nz*nx*sizeof(Float32)
       end
       for it = 1 : nt
           position = (it-1)*sptSize + inds
           seek(fid, position)
           tmp = read(fid, Float32, nz*nx)
           write(fido, tmp)
       end
    end
    close(fido)
    return nothing
end

function NormOfSnapShots(pathin::String, itype="l2")
    fid = open(pathin, "r")
    nz  = read(fid, Int32)
    nx  = read(fid, Int32)
    ext = read(fid, Int32)
    iflag = read(fid, Int32)
    dt    = read(fid, Float32)
    if iflag == 1
       Nz = nz + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    snapSize = sizeof(Float32)*Nz*Nx*10
    nt = round(Int64, ((filesize(fid)-sizeof(Float32)*5)/snapSize))
    if itype == "l2"
       tmp = 0.0
       for it = 1 : nt
           val = convert(Array{Float64}, read(fid, Float32, Nz*Nx*10))
           tmp = tmp + dot(val, val)
       end
       tmp = sqrt(tmp)
    elseif itype == "l1"
       tmp = 0.0
       for it = 1 : nt
           val = convert(Array{Float64}, read(fid, Float32, Nz*Nx*10))
           tmp = tmp + sum(abs(val))
       end
    end
    close(fid)
    return tmp
end

function NormOfSpt(snapShot::SnapShot, itype="l2")
  tmp = 0.0
  if itype == "l2"
     tmp = tmp + dot(snapShot.Vxx, snapShot.Vxx)
     tmp = tmp + dot(snapShot.Vxz, snapShot.Vxz)
     tmp = tmp + dot(snapShot.Vzx, snapShot.Vzx)
     tmp = tmp + dot(snapShot.Vzz, snapShot.Vzz)
     tmp = tmp + dot(snapShot.Txxx, snapShot.Txxx)
     tmp = tmp + dot(snapShot.Txxz, snapShot.Txxz)
     tmp = tmp + dot(snapShot.Tzzx, snapShot.Tzzx)
     tmp = tmp + dot(snapShot.Tzzz, snapShot.Tzzz)
     tmp = tmp + dot(snapShot.Txzx, snapShot.Txzx)
     tmp = tmp + dot(snapShot.Txzz, snapShot.Txzz)
     tmp = sqrt(tmp)
  elseif itype == "l1"
     tmp = tmp + sum(abs(snapShot.Vxx))
     tmp = tmp + sum(abs(snapShot.Vxz))
     tmp = tmp + sum(abs(snapShot.Vzx))
     tmp = tmp + sum(abs(snapShot.Vzz))
     tmp = tmp + sum(abs(snapShot.Txxx))
     tmp = tmp + sum(abs(snapShot.Txxz))
     tmp = tmp + sum(abs(snapShot.Tzzx))
     tmp = tmp + sum(abs(snapShot.Tzzz))
     tmp = tmp + sum(abs(snapShot.Txzx))
     tmp = tmp + sum(abs(snapShot.Txzz))
  end
  return tmp
end

function ScaleWfd!(pathin::String, alpha::Float64)
    alpha = convert(Float32, alpha)
    fid = open(pathin, "r+")
    nz    = read(fid, Int32)
    nx    = read(fid, Int32)
    snapSize = sizeof(Float32)*nz*nx*5
    nt = round(Int64, ((filesize(fid)-sizeof(Float32)*5)/snapSize))
    for it =1: nt
        position = sizeof(Float32)*5 + (it-1)*snapSize
        seek(fid, position)
        tmp = read(fid, Float32, nz*nx*5) .* alpha
        seek(fid, position)
        write(fid, tmp)
    end
    close(fid)
    return nothing
end

function ReverseOrderWfd(pathin::String)
    fid = open(pathin, "r")
    pathin_tmp = join([pathin "_tmp"])
    fid1  = open(pathin_tmp, "w")
    nz    = read(fid, Int32)
    nx    = read(fid, Int32)
    ext   = read(fid, Int32)
    iflag = read(fid, Int32)
    dt    = read(fid, Float32)
    write(fid1, nz, nx, ext, iflag, dt);
    snapSize = sizeof(Float32)*nz*nx*5
    nt = round(Int32, ((filesize(fid)-sizeof(Float32)*5)/snapSize))
    for it = nt : -1 : 1
        position = sizeof(Float32)*5 + (it-1)*snapSize
        seek(fid, position)
        tmp = read(fid, Float32, nz*nx*5)
        write(fid1, tmp)
    end
    close(fid); close(fid1)
    mv(pathin_tmp, pathin, remove_destination=true)
    return nothing
end

function IptWfd(path::String, path1::String)
    fid  = open(path , "r")
    fid1 = open(path1, "r")
    nz =read(fid ,Int32); nx =read(fid ,Int32); ext =read(fid ,Int32); iflag =read(fid ,Int32);
    nz1=read(fid1,Int32); nx1=read(fid1,Int32); ext1=read(fid1,Int32); iflag1=read(fid1,Int32);
    if nz1!=nz || nx1!=nx || ext1!=ext || iflag1!=iflag
       error("size dismatch")
    end
    snapSize = sizeof(Float32)*nz*nx*5
    nt = round(Int32, ((filesize(fid)-sizeof(Float32)*5)/snapSize))
    seek(fid , sizeof(Float32)*5)
    seek(fid1, sizeof(Float32)*5)
    tmp = 0.0
    for it = 1 : nt
        tmp1 = read(fid , Float32, nz*nx*5)
        tmp2 = read(fid1, Float32, nz*nx*5)
        tmp  = tmp + dot(tmp, tmp1)
    end
    close(fid); clsoe(fid1)
    return tmp
end

function partialV(spt::SnapShot, dz::Float64, dx::Float64)
    a1 = 9/8; a2 = -1/24;
    nz = spt.nz; nx = spt.nx; ext = spt.ext; iflag = spt.iflag
    vxx = zeros(nz,nx); vxz = zeros(nz,nx);
    vzx = zeros(nz,nx); vzz = zeros(nz,nx);
    if iflag == 1
       Nz = nz + 2*ext
       upper = ext
    elseif iflag == 2
       Nz = nz +   ext
       upper = 0
    end
    vx = spt.Vxx + spt.Vxz
    vz = spt.Vzx + spt.Vzz
    for iz = 1 : nz
        for ix = 1 : nx
          # compute Dvx/Dx, centered on (i,j)
            ind1 = (ix+ext-3)*Nz + (iz+upper)
            ind2 = (ix+ext-2)*Nz + (iz+upper)
            ind3 = (ix+ext-1)*Nz + (iz+upper)
            ind4 = (ix+ext  )*Nz + (iz+upper)
            vxx[iz,ix] =  a1/dx*(vx[ind3]-vx[ind2]) + a2/dx*(vx[ind4]-vx[ind1])
          # compute Dvx/Dz
            if iflag == 2 && iz == 1
               ind1 = (ix+ext-1)*Nz + iz+upper
               ind2 = (ix+ext-1)*Nz + iz+upper+1
               ind3 = (ix+ext-1)*Nz + iz+upper+2
               vxz[iz,ix] =  a1/dz*(vx[ind2]-vx[ind1]) + a2/dx*(vx[ind3])
            else
               ind1 = (ix+ext-1)*Nz + iz+upper-1
               ind2 = (ix+ext-1)*Nz + iz+upper
               ind3 = (ix+ext-1)*Nz + iz+upper+1
               ind4 = (ix+ext-1)*Nz + iz+upper+2
               vxz[iz,ix] =  a1/dz*(vx[ind3]-vx[ind2]) + a2/dz*(vx[ind4]-vx[ind1])
            end
            # compute Dvz/Dx
               ind1 = (ix+ext-2)*Nz + (iz+upper)
               ind2 = (ix+ext-1)*Nz + (iz+upper)
               ind3 = (ix+ext  )*Nz + (iz+upper)
               ind4 = (ix+ext+1)*Nz + (iz+upper)
               vzx[iz,ix] =  a1/dx*(vz[ind3]-vz[ind2]) + a2/dx*(vz[ind4]-vz[ind1])
          # compute Dvz/Dz
            if iflag == 2 && iz == 1
               ind1 = (ix+ext-1)*Nz + iz+upper
               ind2 = (ix+ext-1)*Nz + iz+upper+1
               vzz[iz,ix] =  a1/dz*vz[ind1] + a2/dz*vz[ind2]
            elseif iflag == 2 && iz == 2
               ind1 = (ix+ext-1)*Nz + iz+upper-1
               ind2 = (ix+ext-1)*Nz + iz+upper
               ind3 = (ix+ext-1)*Nz + iz+upper+1
               vzz[iz,ix] =  a1/dz*(vz[ind2]-vz[ind1]) + a2/dz*vz[ind3]
            else
               ind1 = (ix+ext-1)*Nz + iz+upper-2
               ind2 = (ix+ext-1)*Nz + iz+upper-1
               ind3 = (ix+ext-1)*Nz + iz+upper
               ind4 = (ix+ext-1)*Nz + iz+upper+1
               vzz[iz,ix] =  a1/dz*(vz[ind3]-vz[ind2]) + a2/dz*(vz[ind4]-vz[ind1])
            end
        end
    end
    return vxx, vxz, vzx, vzz
end

function WritePv(path::String, vxx::Array{Float64,2}, vxz::Array{Float64,2}, vzx::Array{Float64,2}, vzz::Array{Float64,2}, nz::Int64, nx::Int64, ext::Int64, iflag::Int64, dt::Float64)
    fid = open(path, "w")
    write(fid, Int32(nz), Int32(nx), Int32(ext), Int32(iflag), Float32(dt))
    write(fid, convert(Array{Float32}, vec(vxx)))
    write(fid, convert(Array{Float32}, vec(vxz)))
    write(fid, convert(Array{Float32}, vec(vzx)))
    write(fid, convert(Array{Float32}, vec(vzz)))
    flush(fid)
    return fid
end

function WritePv(fid::IOStream, vxx::Array{Float64,2}, vxz::Array{Float64,2}, vzx::Array{Float64,2}, vzz::Array{Float64,2})
    write(fid, convert(Array{Float32}, vec(vxx)))
    write(fid, convert(Array{Float32}, vec(vxz)))
    write(fid, convert(Array{Float32}, vec(vzx)))
    write(fid, convert(Array{Float32}, vec(vzz)))
    flush(fid)
    return nothing
end

function ReadPv(path::String, it::Int64)
    fid = open(path, "r")
    nz  = Int64(read(fid,Int32)); nx = Int64(read(fid,Int32));
    ext = Int64(read(fid,Int32)); iflag = Int64(read(fid,Int32));
    dt  = Float64(read(fid,Float32));
    sliceSize = sizeof(Float32) * nz * nx * 4
    position  = sizeof(Int32)*5 + (it-1)*sliceSize
    seek(fid, position)
    vxx = convert(Array{Float64}, read(fid, Float32, nz*nx))
    vxz = convert(Array{Float64}, read(fid, Float32, nz*nx))
    vzx = convert(Array{Float64}, read(fid, Float32, nz*nx))
    vzz = convert(Array{Float64}, read(fid, Float32, nz*nx))
    vxx = reshape(vxx, nz, nx)
    vxz = reshape(vxz, nz, nx)
    vzx = reshape(vzx, nz, nx)
    vzz = reshape(vzz, nz, nx)
    close(fid)
    return vxx, vxz, vzx, vzz
end

function InfoPv(path::String; print_flag=false)
    fid = open(path, "r")
    nz  = Int64(read(fid,Int32)); nx = Int64(read(fid,Int32));
    ext = Int64(read(fid,Int32)); iflag = Int64(read(fid,Int32));
    dt = Float64(read(fid, Float32));
    sliceSize = sizeof(Float32) * nz * nx * 4
    ns = round(Int64, (filesize(fid)-sizeof(Int32)*4)/sliceSize)
    if print_flag
       println("nz   : $nz"   )
       println("nx   : $nx"   )
       println("ext  : $ext"  )
       println("iflag: $iflag")
       println("dt   : $dt"   )
       println("ns   : $ns"   )
    end
    close(fid)
    return nz, nx, ext, iflag, dt, ns
end

function IptSpts(path::String, path1::String)
    (nz, nx, ext, iflag, dt, nt) = InfoSnapShots(path1)
    tmp = 0.0
    for it = 1 : nt
        spt = ReadSnapShot(path , it)
        spt1= ReadSnapShot(path1, it)
        tmp = tmp + dot(spt.Vxx , spt1.Vxx )
        tmp = tmp + dot(spt.Vxz , spt1.Vxz )
        tmp = tmp + dot(spt.Vzx , spt1.Vzx )
        tmp = tmp + dot(spt.Vzz , spt1.Vzz )
        tmp = tmp + dot(spt.Txxx, spt1.Txxx)
        tmp = tmp + dot(spt.Txxz, spt1.Txxz)
        tmp = tmp + dot(spt.Tzzx, spt1.Tzzx)
        tmp = tmp + dot(spt.Tzzz, spt1.Tzzz)
        tmp = tmp + dot(spt.Txzx, spt1.Txzx)
        tmp = tmp + dot(spt.Txzz, spt1.Txzz)
    end
    return tmp
end
