type Shot
    nr   :: Int64
    isz  :: Int64
    isx  :: Int64
    nz   :: Int64
    nx   :: Int64
    ext  :: Int64
    iflag:: Int64
    irz  :: Array{Int64,1}
    irx  :: Array{Int64,1}
    ot   :: Float64
    dt   :: Float64
    nt   :: Int64
    Vx   :: Array{Float64,2}
    Vz   :: Array{Float64,2}
    Txx  :: Array{Float64,2}
    Tzz  :: Array{Float64,2}
    Txz  :: Array{Float64,2}
end

function InitShot(isz::Int64, isx::Int64, nz::Int64, nx::Int64, ext::Int64, iflag::Int64, irz::Array{Int64,1}, irx::Array{Int64,1}, ot::Float64, dt::Float64, nt::Int64)
    nr = length(irz)
    shot = Shot(nr, isz, isx, nz, nx, ext, iflag, irz, irx, ot, dt, nt, zeros(nt,nr), zeros(nt,nr), zeros(nt,nr), zeros(nt,nr), zeros(nt,nr))
    return shot
end

type ShotV
    nr   :: Int64
    isz  :: Int64
    isx  :: Int64
    nz   :: Int64
    nx   :: Int64
    ext  :: Int64
    iflag:: Int64
    irz  :: Array{Int64,1}
    irx  :: Array{Int64,1}
    ot   :: Float64
    dt   :: Float64
    nt   :: Int64
    Vx   :: Array{Float64,2}
    Vz   :: Array{Float64,2}
end

function InitShotV(isz::Int64, isx::Int64, nz::Int64, nx::Int64, ext::Int64, iflag::Int64, irz::Array{Int64,1}, irx::Array{Int64,1}, ot::Float64, dt::Float64, nt::Int64)
    nr = length(irz)
    shotv = ShotV(nr, isz, isx, nz, nx, ext, iflag, irz, irx, ot, dt, nt, zeros(nt,nr), zeros(nt,nr))
    return shotv
end

function WriteShot(pathout::String, shot::Shot)
    fid = open(pathout, "w")
    nr  = length(shot.irz)
    write(fid, Int32(nr),      Int32(shot.isz), Int32(shot.isx))
    write(fid, Int32(shot.nz), Int32(shot.nx) , Int32(shot.ext), Int32(shot.iflag))
    write(fid, convert(Array{Int32},shot.irz) , convert(Array{Int32},shot.irx))
    write(fid, Float32(shot.ot), Float32(shot.dt), Int32(shot.nt))
    write(fid, convert(Array{Float32},vec(shot.Vx)), convert(Array{Float32},vec(shot.Vz)), convert(Array{Float32},vec(shot.Txx)), convert(Array{Float32},vec(shot.Tzz)), convert(Array{Float32},vec(shot.Txz)))
    close(fid)
    return nothing
end

function WriteShotV(path::String, shot::ShotV)
    fid = open(path, "w")
    nr  = length(shot.irz)
    write(fid, Int32(nr),      Int32(shot.isz), Int32(shot.isx))
    write(fid, Int32(shot.nz), Int32(shot.nx) , Int32(shot.ext), Int32(shot.iflag))
    write(fid, convert(Array{Int32},shot.irz) , convert(Array{Int32},shot.irx))
    write(fid, Float32(shot.ot), Float32(shot.dt), Int32(shot.nt))
    write(fid, convert(Array{Float32},vec(shot.Vx)), convert(Array{Float32},vec(shot.Vz)))
    close(fid)
    return nothing
end

function ReadShot(pathin::String)
    fid = open(pathin, "r")
    nr  = Int64(read(fid,Int32)); isz = Int64(read(fid,Int32)); isx = Int64(read(fid,Int32));
    nz  = Int64(read(fid,Int32)); nx  = Int64(read(fid,Int32)); ext = Int64(read(fid,Int32)); iflag = Int64(read(fid,Int32));
    irz = convert(Array{Int64},read(fid,Int32,nr)); irx = convert(Array{Float64},read(fid,Int32,nr));
    ot  = Float64(read(fid,Float32))  ; dt  = Float64(read(fid,Float32))  ; nt = Int64(read(fid,Int32));
    Vx  = convert(Array{Float64}, reshape(read(fid,Float32,nt*nr),nt,nr))
    Vz  = convert(Array{Float64}, reshape(read(fid,Float32,nt*nr),nt,nr))
    Txx = convert(Array{Float64}, reshape(read(fid,Float32,nt*nr),nt,nr))
    Tzz = convert(Array{Float64}, reshape(read(fid,Float32,nt*nr),nt,nr))
    Vxz = convert(Array{Float64}, reshape(read(fid,Float32,nt*nr),nt,nr))
    close(fid)
    shot = Shot(nr, isz, isx, nz, nx, ext, iflag, irz, irx, ot, dt, nt, Vx, Vz, Txx, Tzz, Txz)
    return shot
end

function ReadShotV(pathin::String)
    fid = open(pathin, "r")
    nr  = Int64(read(fid,Int32)); isz = Int64(read(fid,Int32)); isx = Int64(read(fid,Int32));
    nz  = Int64(read(fid,Int32)); nx  = Int64(read(fid,Int32)); ext = Int64(read(fid,Int32)); iflag = Int64(read(fid,Int32));
    irz = convert(Array{Int64},read(fid,Int32,nr)); irx = convert(Array{Int64},read(fid,Int32,nr));
    ot  = Float64(read(fid,Float32))  ; dt  = Float64(read(fid,Float32))  ; nt = Int64(read(fid,Int32));
    Vx  = convert(Array{Float64}, reshape(read(fid,Float32,nt*nr),nt,nr))
    Vz  = convert(Array{Float64}, reshape(read(fid,Float32,nt*nr),nt,nr))
    close(fid)
    shot = ShotV(nr, isz, isx, nz, nx, ext, iflag, irz, irx, ot, dt, nt, Vx, Vz)
    return shot
end

function InfoShotV(pathin::String)
    fid = open(pathin, "r")
    nr  = Int64(read(fid,Int32)); isz = Int64(read(fid,Int32)); isx = Int64(read(fid,Int32));
    nz  = Int64(read(fid,Int32)); nx  = Int64(read(fid,Int32)); ext = Int64(read(fid,Int32)); iflag = Int64(read(fid,Int32));
    irz = convert(Array{Int64},read(fid,Int32,nr)); irx = convert(Array{Int64},read(fid,Int32,nr));
    ot  = Float64(read(fid,Float32))  ; dt  = Float64(read(fid,Float32))  ; nt = Int64(read(fid,Int32));
    return nr, nz, nx, ext, iflag, irz, irx, dt, nt
end

function Spt2Shot!(shot::Shot, spt::SnapShot)
    it = spt.it;
    nz=shot.nz; nx=shot.nx; ext=shot.ext; iflag=shot.iflag;
    if iflag == 1
       Nz = nz + 2*ext
       upper = ext
    elseif iflag == 2
       Nz = nz +   ext
       upper = 0
    end
    nr = length(shot.irz)
    for ir = 1: nr
        iz = shot.irz[ir]
        ix = shot.irx[ir]
        ind = (ext+ix-1)*Nz + (iz+upper)
        shot.Vx[it,ir]  = spt.Vxx[ind] + spt.Vxz[ind]
        shot.Vz[it,ir]  = spt.Vzx[ind] + spt.Vzz[ind]
        shot.Txx[it,ir] = spt.Txxx[ind] + spt.Txxz[ind]
        shot.Tzz[it,ir] = spt.Tzzx[ind] + spt.Tzzz[ind]
        shot.Txz[it,ir] = spt.Txzx[ind] + spt.Txzz[ind]
    end
    return nothing
end

function Spt2ShotV!(shotv::ShotV, spt::SnapShot)
    it = spt.it;
    nz=shotv.nz; nx=shotv.nx; ext=shotv.ext; iflag=shotv.iflag;
    if iflag == 1
       Nz = nz + 2*ext
       upper = ext
    elseif iflag == 2
       Nz = nz +   ext
       upper = 0
    end
    nr = length(shotv.irz)
    for ir = 1: nr
        iz = shotv.irz[ir]
        ix = shotv.irx[ir]
        ind = (ext+ix-1)*Nz + (iz+upper)
        shotv.Vx[it,ir]  = spt.Vxx[ind] + spt.Vxz[ind]
        shotv.Vz[it,ir]  = spt.Vzx[ind] + spt.Vzz[ind]
    end
    return nothing
end

function AddShot2Spt!(spt::SnapShot, shot::Shot)
    it = spt.it
    if spt.iflag == 1
       Nz = spt.nz + 2*spt.ext
       upper = spt.ext
    elseif spt.iflag == 2
       Nz = spt.nz +   spt.ext
       upper = 0
    end
    nr = length(shot.irz)
    for ir = 1: nr
        iz = shot.irz[ir]
        ix = shot.irx[ir]
        ind = (spt.ext+ix-1)*Nz + (iz+upper)
        spt.Vxx[ind] = spt.Vxx[ind] + shot.Vx[it,ir]
        spt.Vxz[ind] = spt.Vxz[ind] + shot.Vx[it,ir]
        spt.Vzx[ind] = spt.Vzx[ind] + shot.Vz[it,ir]
        spt.Vzz[ind] = spt.Vzz[ind] + shot.Vz[it,ir]
        spt.Txxx[ind] = spt.Txxx[ind] + shot.Txx[it,ir]
        spt.Txxz[ind] = spt.Txxz[ind] + shot.Txx[it,ir]
        spt.Tzzx[ind] = spt.Tzzx[ind] + shot.Tzz[it,ir]
        spt.Tzzz[ind] = spt.Tzzz[ind] + shot.Tzz[it,ir]
        spt.Txzx[ind] = spt.Txzx[ind] + shot.Txz[it,ir]
        spt.Txzz[ind] = spt.Txzz[ind] + shot.Txz[it,ir]
    end
    return nothing
end

function AddShotV2Spt!(spt::SnapShot, shotv::ShotV)
    it = spt.it
    if spt.iflag == 1
       Nz = spt.nz + 2*spt.ext
       upper = spt.ext
    elseif spt.iflag == 2
       Nz = spt.nz +   spt.ext
       upper = 0
    end
    nr = length(shotv.irz)
    for ir = 1: nr
        iz = shotv.irz[ir]
        ix = shotv.irx[ir]
        ind = (spt.ext+ix-1)*Nz + (iz+upper)
        spt.Vxx[ind] = spt.Vxx[ind] + shotv.Vx[it,ir]
        spt.Vxz[ind] = spt.Vxz[ind] + shotv.Vx[it,ir]
        spt.Vzx[ind] = spt.Vzx[ind] + shotv.Vz[it,ir]
        spt.Vzz[ind] = spt.Vzz[ind] + shotv.Vz[it,ir]
    end
    return nothing
end

function Dif2Records(records1, records2) #save result of records2 - records1 to records1
  if records1[1].nz != records2[1].nz || records1[1].nx != records2[1].nx || records1[1].ext != records2[1].ext
     error("size does not match")
  end
  if records1[1].nt != records2[1].nt
     error("time length does not match")
  end
  nrec = length(records1)
  for irec = 1: nrec
      records1[irec].Vx =  records2[irec].Vx -  records1[irec].Vx
      records1[irec].Vz =  records2[irec].Vz -  records1[irec].Vz
      records1[irec].Txx = records2[irec].Txx - records1[irec].Txx
      records1[irec].Tzz = records2[irec].Tzz - records1[irec].Tzz
      records1[irec].Txz = records2[irec].Txz - records1[irec].Txz
  end
  return records1
end

function Dif2RecordsNew(records1, records2) #save result of records1 - records2 to records2
  if records1[1].nz != records2[1].nz || records1[1].nx != records2[1].nx || records1[1].ext != records2[1].ext
     error("size does not match")
  end
  if records1[1].nt != records2[1].nt
     error("time length does not match")
  end
  nrec = length(records1)
  for irec = 1: nrec
      records2[irec].Vx =  records1[irec].Vx -  records2[irec].Vx
      records2[irec].Vz =  records1[irec].Vz -  records2[irec].Vz
      records2[irec].Txx = records1[irec].Txx - records2[irec].Txx
      records2[irec].Tzz = records1[irec].Tzz - records2[irec].Tzz
      records2[irec].Txz = records1[irec].Txz - records2[irec].Txz
  end
  return records2
end


function NormOfShot(shpt::Shot; itype="l2")
  tmp = 0.0
  nr  = length(shot.irz)
  if itype == "l2"
     tmp = tmp + dot(vec(shot.Vx), vec(shot.Vx )) + dot(vec(shot.Vz ),vec(shot.Vz ))
               + dot(vec(shot.Txx),vec(shot.Txx)) + dot(vec(shot.Tzz),vec(shot.Tzz)) + dot(vec(shot.Txz),vec(shot.Txz))
     tmp = sqrt(tmp)
  elseif itype == "l1"
     tmp = tmp + sum(abs(shot.Vx )) + sum(abs(shot.Vz ))
               + sum(abs(shot.Txx)) + sum(abs(shot.Tzz)) + sum(abs(shot.Txz))
  end
  return tmp
end

function imaging(dm::Array{Float64,2}, du::Array{Float64,2}, spt::SnapShot, path::String)
    nz = spt.nz; nx = spt.nx; ext = spt.ext; iflag = spt.iflag;
    it = spt.it
    if spt.iflag == 1
       Nz = spt.nz + 2*spt.ext
       upper = spt.ext
    elseif spt.iflag == 2
       Nz = spt.nz +   spt.ext
       upper = 0
    end
    Nx = spt.nx + 2*spt.ext
    (vxx, vxz, vzx, vzz) = ReadPv(path, it)
    Txx = reshape(spt.Txxx+spt.Txxz, Nz, Nx) * 1/2
    Tzz = reshape(spt.Tzzx+spt.Tzzz, Nz, Nx) * 1/2
    Txz = reshape(spt.Txzx+spt.Txzz, Nz, Nx) * 1/2
    tmp = Txx[upper+1:upper+nz, ext+1:ext+nx] + Tzz[upper+1:upper+nz, ext+1:ext+nx]
    dm  = dm + (vxx+vzz) .* tmp
    du  = du + 2*vxx.*Txx[upper+1:upper+nz, ext+1:ext+nx]
    du  = du + 2*vzz.*Tzz[upper+1:upper+nz, ext+1:ext+nx]
    du  = du + (vxz+vzx) .* Txz[upper+1:upper+nz, ext+1:ext+nx]
    return dm, du
end

function SptsPv2Born(path_spt, path_pv)
    (nz, nx, ext, iflag, dt, nt) = InfoPv(path_pv)
    if iflag == 1
       Nz = nz + 2*ext
       upper = ext
    elseif iflag == 2
       Nz = nz +   ext
       upper = 0
    end
    Nx = nx + 2*ext
    dm = zeros(nz, nx); du = zeros(nz, nx);
    for it = 1 : nt
        (vxx, vxz, vzx, vzz) = ReadPv(path_pv, it)
        spt = ReadSnapShot(path_spt, it)
        Txx = reshape(spt.Txxx+spt.Txxz, Nz, Nx) * 1/2
        Tzz = reshape(spt.Tzzx+spt.Tzzz, Nz, Nx) * 1/2
        Txz = reshape(spt.Txzx+spt.Txzz, Nz, Nx) * 1/2
        tmp = Txx[upper+1:upper+nz, ext+1:ext+nx] + Tzz[upper+1:upper+nz, ext+1:ext+nx]
        dm  = dm + (vxx+vzz) .* tmp
        du  = du + 2*vxx.*Txx[upper+1:upper+nz, ext+1:ext+nx]
        du  = du + 2*vzz.*Tzz[upper+1:upper+nz, ext+1:ext+nx]
        du  = du + (vxz+vzx) .* Txz[upper+1:upper+nz, ext+1:ext+nx]
    end
    return dm, du
end

function removeDirect(shotv::ShotV, v::Float64, src::Source, dx::Float64, dt::Float64)
    nr = length(shotv.irz)
    nt = size(shotv.Vx, 1 )
    isz = src.isz; isx = src.isx;
    lsrc= length(src.Txx)
    damp= hanning(21)[2:11]
    if shotv.irz[1] != isz
       error("source and receivers are not on same level")
    end
    for ir = 1 : nr
        dp  = ones(nt)
        l  = abs((shotv.irx[ir]-isx)*dx)
        it = ceil(Int64,l/(v*dt))+1+lsrc
        dp[1:it] = zeros(it)
        dp[it+1:it+10] = damp
        shotv.Vx[:,ir] = shotv.Vx[:,ir] .* dp
        shotv.Vz[:,ir] = shotv.Vz[:,ir] .* dp
    end
    return nothing
end

function binWfd2shotv(irz::Array{Int64,1}, irx::Array{Int64,1}, path::String)
    fid = open(path, "r")
    nr  = length(irz)
    nz = Int64(read(fid, Int32))
    nx = Int64(read(fid, Int32))
    ext= Int64(read(fid, Int32))
    iflag= Int64(read(fid, Int32))
    dt = Float64(read(fid, Float32))
    lwfd = sizeof(Float32)*nz*nx*5
    pre  = sizeof(Float32)*5
    nt = round(Int64, filesize(fid)/lwfd)
    shotv = InitShotV(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
    for it = 1 : nt
        pos = pre + (it-1)*lwfd
        seek(fid, pos)
        vx = reshape(read(fid, Float32, nz*nx), nz, nx)
        vz = reshape(read(fid, Float32, nz*nx), nz, nx)
        for ir = 1 : nr
            shotv.Vx[it,ir] = vx[irz[ir], irx[ir]]
            shotv.Vz[it,ir] = vz[irz[ir], irx[ir]]
        end
    end
    return shotv
end
