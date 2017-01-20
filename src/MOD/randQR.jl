type Psrc
     k  :: Int64
     ot :: Float64
     dt :: Float64
     nt :: Int64
     sc :: Array{Float64, 2}
end

type ModSpt
     k  :: Int64
     dt :: Float64
     it :: Int64
     wc :: Array{Float64,1}
end

function InitModSpt(k::Int64, dt::Float64, it::Int64)
    spt = ModSpt(k, dt, it, zeros(k))
    return spt
end

function ReadModSpt(path::String, it::Int64)
    fid = open(path, "r")
    k = Int64(read(fid, Int32))
    dt= Float64(read(fid, Float32))
    nt= floor(Int64, (filesize(fid)-sizeof(Int32)*2)/(sizeof(Float32)*k))
    pos = sizeof(Float32)*2 + sizeof(Float32)*k;
    seek(fid, pos)
    wc = convert(Array{Float64}, read(fid, Float32, k))
    return ModSpt(k, dt, it, wc)
end

function WriteModSpt(path::String, spt::ModSpt)
    fid = open(path, "w")
    write(fid, Int32(spt.k), Float32(spt.dt))
    write(fid, convert(Array{Float32}, spt.wc))
    flush(fid)
    return fid
end

function WriteModSpt(fid::IOStream, spt::ModSpt)
    write(fid, convert(Array{Float32}, spt.wc))
    flush(fid)
    return nothing
end

function InfoModSpts(path::String)
    fid = open(path, "r")
    k = Int64(read(fid, Int32))
    dt= Float64(read(fid, Float32))
    nt= floor(Int64, (filesize(fid)-sizeof(Int32)*2)/(sizeof(Float32)*k))
    return k, dt, nt
end

function MultiStepForward_MOR(pathc::String, pathFr::String, spj::Psrc; tmax=1.0)
    Fr = ReadRDFD(pathFr); k = size(Fr, 1);
    dt = spj.dt
    tl = spj.ot; tu = tl + (spj.nt-1)*dt
    nt = floor(Int64, tmax/dt) + 1
    spt1 = InitModSpt(k, dt, 1)
    spt2 = InitModSpt(k, dt, 2)
    AddPsrc!(spt1, spj)
    fid = WriteModSpt(pathc, spt1)
    for it = 2 : nt
        spt2.wc=Fr*spt1.wc; spt2.it=spt1.it+1;
        if tl <= (it-1)*dt <= tu
           AddPsrc!(spt2, spj)
        end
        spt1.wc[:] = spt2.wc[:]; spt1.it = spt2.it;
        WriteModSpt(fid, spt1)
    end
    close(fid)
end

function ModSpts2Spts(path::String, pathe::String, pathc::String; rk=0)
    (nz, nx, ext, iflag, dt, nt) = InfoSnapShots(pathe)
    if rk == 0
       rk = nt
    end
    Q = formBaseMatrix(pathe, rk=rk)
    Q = convert(Array{Float32}, Q)
    fid = open(path, "w")
    write(fid, Int32(nz), Int32(nx), Int32(ext), Int32(iflag), Float32(dt))
    (k, dt, nt) = InfoModSpts(pathc)
    fidc = open(pathc, "r"); seek(fidc, sizeof(Int32)*2)
    for it = 1 : nt
        wc = read(fidc, Float32, k)
        spt = Q * wc
        write(fid, spt)
    end
    close(fidc); close(fid);
    return nothing
end

function SrcProj(src::Source, pathe::String; rk=0)
    Qt = formBaseMatrix(pathe, rk=rk); Qt = Qt';
    k  = size(Qt, 1)
    (tl, tu) = SrcRange(src)
    dt = src.dt; nt = round(Int64, (tu-tl)/dt) + 1;
    spj= Psrc(k, tl, dt, nt, zeros(Float64, k, nt))
    nz = src.nz; nx = src.nx; ext= src.ext; iflag = src.iflag;
    if iflag == 1
       Nz = nz + 2*ext
       upper = ext
    elseif iflag == 2
       Nz = nz +   ext
       upper = 0
    end
    Nx = nx + 2*ext
    iz = src.isz + upper
    ix = src.isx + ext
    for it = 1 : nt
        if src.Vx_flag
           ind= (ix-1)*Nz+iz
           spj.sc[:,it] = spj.sc[:,it] + Qt[:,ind]*src.Vx[it]/2
           ind= (ix-1)*Nz+iz + Nz*Nx*1
           spj.sc[:,it] = spj.sc[:,it] + Qt[:,ind]*src.Vx[it]/2
        end
        if src.Vz_flag
           ind= (ix-1)*Nz+iz + Nz*Nx*2
           spj.sc[:,it] = spj.sc[:,it] + Qt[:,ind]*src.Vz[it]/2
           ind= (ix-1)*Nz+iz + Nz*Nx*3
           spj.sc[:,it] = spj.sc[:,it] + Qt[:,ind]*src.Vz[it]/2
        end
        if src.Txx_flag
           ind= (ix-1)*Nz+iz + Nz*Nx*4
           spj.sc[:,it] = spj.sc[:,it] + Qt[:,ind]*src.Txx[it]/2
           ind= (ix-1)*Nz+iz + Nz*Nx*5
           spj.sc[:,it] = spj.sc[:,it] + Qt[:,ind]*src.Txx[it]/2
        end
        if src.Tzz_flag
           ind= (ix-1)*Nz+iz + Nz*Nx*6
           spj.sc[:,it] = spj.sc[:,it] + Qt[:,ind]*src.Tzz[it]/2
           ind= (ix-1)*Nz+iz + Nz*Nx*7
           spj.sc[:,it] = spj.sc[:,it] + Qt[:,ind]*src.Tzz[it]/2
        end
        if src.Txz_flag
           ind= (ix-1)*Nz+iz + Nz*Nx*8
           spj.sc[:,it] = spj.sc[:,it] + Qt[:,ind]*src.Txz[it]/2
           ind= (ix-1)*Nz+iz + Nz*Nx*9
           spj.sc[:,it] = spj.sc[:,it] + Qt[:,ind]*src.Txz[it]/2
        end
    end
    return spj
end

function AddPsrc!(spt::ModSpt, spj::Psrc)
    dt = spt.dt
    it = spt.it
    if spj.ot <= (it-1)*dt <= spj.ot + (spj.nt-1)*dt
       indt = it - round(Int64, spj.ot/dt)
       spt.wc[:] = spt.wc[:] + spj.sc[:,indt]
    end
    return nothing
end

function lowRankFd(path_Fr::String, fidMtx::FidMtx, pathe::String; rk=0)
    (nz, nx, ext, iflag, dt, nt) = InfoSnapShots(pathe)
    path_tmp = join([pathe "_OneStepForwardBase"])
    spt1 = InitSnapShot(nz, nx, ext, iflag, dt, 1)
    spt2 = InitSnapShot(nz, nx, ext, iflag, dt, 2)
    spt1 = ReadSnapShot(pathe, 1)
    OneStepForward!(spt2, spt1, fidMtx)
    fid  = WriteSnapShot(path_tmp, spt2)
    if rk == 0
       rk = nt
    end
    nt = rk < nt ? rk : nt
    for it = 2 : nt
        spt1 = ReadSnapShot(pathe, it)
        OneStepForward!(spt2, spt1, fidMtx)
        WriteSnapShot(fid, spt2)
    end
    close(fid)
    Q = formBaseMatrix(pathe, rk=rk)
    Fr = QFQ(Q, path_tmp); rm(path_tmp);
    fid = open(path_Fr, "w"); k = size(Fr, 1);
    write(fid, k); write(fid, vec(Fr)); close(fid);
    return nothing
end

function QFQ(Q::Array{Float64,2}, path::String)
    (m, rk) = size(Q)
    Fr = zeros(Float64, rk, rk)
    fid = open(path, "r")
    for j = 1 : rk
        pos = sizeof(Float32)*5 + sizeof(Float32)*m*(j-1)
        seek(fid, pos)
        tmp = convert(Array{Float64}, read(fid, Float32, m))
        for i = 1 : rk
            Fr[i,j] = dot(Q[:,i], tmp)
        end
    end
    return Fr
end

function ReadRDFD(path::String)
    fid = open(path, "r")
    k = read(fid, Int64)
    Fr = reshape(read(fid, Float64, k*k), k, k)
    return Fr
end

function QRbase(path::String, path1::String; rk=200, style="mem")
    randMix(path, path1, rk=rk, style=style)
    Q = formStateMatrix(path); rm(path);
    lambda = QR(Q)
    (nz, nx, ext, iflag, dt, nt) = InfoSnapShots(path1)
    WriteBase(path, Q, nz, nx, ext, iflag, dt)
    return lambda
end

function QR(A; delta=1e-10)
     nt = size(A, 2)
     lambda = zeros(typeof(A[1]), nt)
     for j = 1 : nt
         q = A[:,j]
         for k = 1 : j-1
             a = A[:, k]
             q = q - dot(a, q)*a
         end
         lambda[j] = norm(q)
         q = q / (lambda[j]+delta)
         A[:, j] = q
     end
     return lambda
end

function randMix(path1::String, path::String; rk=200, style="mem")
    (nz, nx, ext, iflag, dt, nt) = InfoSnapShots(path)
    fid = open(path1, "w")
    write(fid, Int32(nz), Int32(nx), Int32(ext), Int32(iflag), Float32(dt))
    if style == "mem"
       d = formStateMatrix(path)
       (m, n) = size(d)
       for it = 1 : rk
           r = convert(Array{Float32}, randn(n))
           tmp = d * r
           write(fid, tmp)
       end
    end
    if style == "hard"
       if iflag == 1
          Nz = nz + 2*ext
       elseif iflag == 2
          Nz = nz +   ext
       end
       Nx = nx + 2*ext
       len = Nz*Nx*10
       S = Array(Float32, len, rk)
       for l = 1 : len
           tmp = zeros(nt)
           for it = 1 : nt
               pos = ((it-1)*len+l-1)*sizeof(Float32)
               seek(fid, pos)
               tmp[it] = read(fid, Float32)
           end
           for k = 1 : rk
               v = convert(Array{Float32}, randn(nt))
               S[l, k] = dot(tmp, v)
           end
       end
       write(fid, S[:])
    end
    close(fid)
end

function formBaseMatrix(pathe::String; rk=0)
    (nz, nx, ext, iflag, dt, nt) = InfoSnapShots(pathe)
    if iflag == 1
       Nz = nz + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    l = Nz * Nx * 10
    if rk == 0
       rk = nt
    end
    if rk > nt
       error("rk can not be larger than nt")
    end
    Q = Array(Float64, l, rk)
    fid = open(pathe, "r")
    sptSize = sizeof(Float32)*Nz*Nx*10
    for it = 1 : rk
        pos = sizeof(Float32)*5 + (it-1)*sptSize
        seek(fid, pos)
        Q[:, it] = convert(Array{Float64}, read(fid,Float32,Nz*Nx*10))
    end
    return Q
end

function formStateMatrix(path::String)
    (nz, nx, ext, iflag, dt, nt) = InfoSnapShots(path)
    if iflag == 1
       Nz = nz + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    l = Nz * Nx * 10
    S = Array(Float32, l, nt)
    fid = open(path, "r")
    sptSize = sizeof(Float32)*Nz*Nx*10
    for it = 1 : nt
        pos = sizeof(Float32)*5 + (it-1)*sptSize
        seek(fid, pos)
        S[:, it] = read(fid, Float32, Nz*Nx*10)
    end
    return S
end

function WriteBase(path::String, Q, nz, nx, ext, iflag, dt)
    fid = open(path, "w")
    write(fid, Int32(nz), Int32(nx), Int32(ext), Int32(iflag), Float32(dt))
    write(fid, convert(Array{Float32}, vec(Q)))
    close(fid)
end
