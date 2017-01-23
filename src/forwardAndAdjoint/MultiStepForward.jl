# return records (vx, vz), inject single source
function MultiStepForward(irz::Array{Int64,1}, irx::Array{Int64,1}, src::Source, fidMtx::FidMtx; tmax=0.5)
    nz  =  src.nz;  nx = src.nx    ;
    ext = src.ext;  iflag=src.iflag; dt = src.dt;
    stl = src.ot ;  stu = src.ot+(src.nt-1)*dt;
    nt  = round(Int64, tmax/dt)+1
    spt1 = InitSnapShot(nz, nx, ext, iflag, dt, 1)
    spt2 = InitSnapShot(nz, nx, ext, iflag, dt, 2)
    AddSource!(spt1, src)
    shotv = InitShotV(src.isz, src.isx, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
    Spt2ShotV!(shotv, spt1)
    tmp = zeros(length(spt1.Vxx))
    tmp1= zeros(tmp)
    for it = 2 : nt
        OneStepForward!(spt2, spt1, fidMtx, tmp, tmp1)
        if stl <= (it-1)*dt <= stu
           AddSource!(spt2, src)
        end
        CopySnapShot!(spt1, spt2)
        Spt2ShotV!(shotv, spt1)
    end
    return shotv
end

# return records (vx, vz), inject multiple sources
function MultiStepForward(irz::Array{Int64,1}, irx::Array{Int64,1}, srcs::Array{Source,1}, fidMtx::FidMtx; tmax=0.5)
    nz  = srcs[1].nz ; nx = srcs[1].nx    ;
    ext = srcs[1].ext; iflag=srcs[1].iflag; dt = srcs[1].dt;
    (stl, stu) = SourcesTimeRange(srcs)
    nt  = round(Int64, tmax/dt)+1
    spt1 = InitSnapShot(nz, nx, ext, iflag, dt, 1)
    spt2 = InitSnapShot(nz, nx, ext, iflag, dt, 2)
    AddMultiSources!(spt1, srcs)
    shotv = InitShotV(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
    Spt2ShotV!(shotv, spt1)
    tmp = zeros(length(spt1.Vxx))
    tmp1= zeros(tmp)
    for it = 2 : nt
        OneStepForward!(spt2, spt1, fidMtx, tmp, tmp1)
        if stl <= (it-1)*dt <= stu
           AddMultiSources!(spt2, srcs)
        end
        CopySnapShot!(spt1, spt2)
        Spt2ShotV!(shotv, spt1)
    end
    return shotv
end

# write wave field(default, 5 components) or snapshot(10 components), inject single source
function MultiStepForward(path::String, src::Source, fidMtx::FidMtx; tmax=0.5, otype="wfd")
    nz  =  src.nz;  nx = src.nx    ;
    ext = src.ext;  iflag=src.iflag; dt  = src.dt;
    stl = src.ot ;  stu = src.ot+(src.nt-1)*dt;
    nt   = round(Int64, tmax/dt)+1
    spt1 = InitSnapShot(nz, nx, ext, iflag, dt, 1)
    spt2 = InitSnapShot(nz, nx, ext, iflag, dt, 2)
    AddSource!(spt1, src)
    if otype == "wfd"
       fid = WriteWfd(path, spt1)
    elseif otype == "spt"
       fid = WriteSnapShot(path, spt1)
    end
    tmp = zeros(length(spt1.Vxx))
    tmp1= zeros(tmp)
    for it = 2 : nt
        OneStepForward!(spt2, spt1, fidMtx, tmp, tmp1)
        if stl <= (it-1)*dt <= stu
           AddSource!(spt2, src)
        end
        CopySnapShot!(spt1, spt2)
        if otype == "wfd"
           WriteWfd(fid, spt1)
        elseif otype == "spt"
           WriteSnapShot(fid, spt1)
        end
    end
    close(fid)
    return nothing
end

# write full wavefield, inject multiple souces
function MultiStepForward(path::String, srcs::Array{Source,1}, fidMtx::FidMtx; tmax=0.5, otype="wfd")
    nz  = srcs[1].nz ; nx   = srcs[1].nx;
    ext = srcs[1].ext; iflag= srcs[1].iflag; dt = srcs[1].dt;
    (stl, stu) = SourcesTimeRange(srcs)
    nt   = round(Int64, tmax/dt)+1
    spt1 = InitSnapShot(nz, nx, ext, iflag, dt, 1)
    spt2 = InitSnapShot(nz, nx, ext, iflag, dt, 2)
    AddMultiSources!(spt1, srcs)
    if otype == "wfd"
       fid = WriteWfd(path, spt1)
    elseif otype == "spt"
       fid = WriteSnapShot(path, spt1)
    end
    tmp = zeros(length(spt1.Vxx))
    tmp1= zeros(tmp)
    for it = 2 : nt
        OneStepForward!(spt2, spt1, fidMtx, tmp, tmp1)
        if stl <= (it-1)*dt <= stu
           AddMultiSources!(spt2, srcs)
        end
        CopySnapShot!(spt1, spt2)
        if otype == "wfd"
           WriteWfd(fid, spt1)
        elseif otype == "spt"
           WriteSnapShot(fid, spt1)
        end
    end
    close(fid)
    return nothing
end

# write (vxx, vxz, vzx, vzz); inject single source
function MultiStepForward(path::String, src::Source, fidMtx::FidMtx, dz::Float64, dx::Float64; tmax=0.5)
    nz  =  src.nz;  nx = src.nx    ;
    ext = src.ext;  iflag=src.iflag; dt  = src.dt;
    stl = src.ot ;  stu = src.ot+(src.nt-1)*dt;
    nt  = round(Int64, tmax/dt)+1
    spt1= InitSnapShot(nz, nx, ext, iflag, dt, 1)
    spt2= InitSnapShot(nz, nx, ext, iflag, dt, 2)
    AddSource!(spt1, src)
    vxx=zeros(nz,nx); vxz=zeros(nz,nx);
    vzx=zeros(nz,nx); vzz=zeros(nz,nx);
    tmp = zeros(length(spt1.Vxx))
    tmp1= zeros(tmp)
    partialV!(vxx, vxz, vzx, vzz, spt1, dz, dx, tmp, tmp1);
    fid = WritePv(path, vxx, vxz, vzx, vzz, nz, nx, ext, iflag, dt);
    for it = 2 : nt
        OneStepForward!(spt2, spt1, fidMtx, tmp, tmp1)
        if stl <= (it-1)*dt <= stu
           AddSource!(spt2, src)
        end
        CopySnapShot!(spt1, spt2)
        partialV!(vxx, vxz, vzx, vzz, spt1, dz, dx, tmp, tmp1); WritePv(fid, vxx, vxz, vzx, vzz);
    end
    close(fid)
    return nothing
end

# write (vxx, vxz, vzx, vzz); inject multiple souces
function MultiStepForward(path::String, srcs::Array{Source,1}, fidMtx::FidMtx, dz::Float64, dx::Float64; tmax=0.5)
    nz  = srcs[1].nz ; nx   = srcs[1].nx   ;
    ext = srcs[1].ext; iflag= srcs[1].iflag; dt = srcs[1].dt;
    (stl, stu) = SourcesTimeRange(srcs)
    nt  = round(Int64, tmax/dt)+1
    spt1= InitSnapShot(nz, nx, ext, iflag, dt, 1)
    spt2= InitSnapShot(nz, nx, ext, iflag, dt, 2)
    AddMultiSources!(spt1, srcs)
    tmp = zeros(length(spt1.Vxx))
    tmp1= zeros(tmp)
    partialV!(vxx, vxz, vzx, vzz, spt1, dz, dx, tmp, tmp1);
    fid = WritePv(path, vxx, vxz, vzx, vzz, nz, nx, ext, iflag, dt);
    for it = 2 : nt
        OneStepForward!(spt2, spt1, fidMtx, tmp, tmp1)
        if stl <= (it-1)*dt <= stu
           AddMultiSources!(spt2, srcs)
        end
        CopySnapShot!(spt1, spt2)
        partialV!(vxx, vxz, vzx, vzz, spt1, dz, dx, tmp, tmp1); WritePv(fid, vxx, vxz, vzx, vzz);
    end
    close(fid)
    return nothing
end

# output shotv(Vx, Vz), Born approximation
function MultiStepForward(irz::Array{Int64,1}, irx::Array{Int64,1}, path::String, dm::Array{Float64,2}, du::Array{Float64,2}, fidMtx::FidMtx; tmax=0.5)
    (nz, nx, ext, iflag, dt, ns) = InfoPv(path);
    nt  = round(Int64, tmax/dt)+1
    spt1= InitSnapShot(nz, nx, ext, iflag, dt, 1)
    spt2= InitSnapShot(nz, nx, ext, iflag, dt, 2)
    AddSourceBorn!(spt1, path, dm, du)
    shotv= InitShotV(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
    Spt2ShotV!(shotv, spt1)
    tmp = zeros(length(spt1.Vxx))
    tmp1= zeros(tmp)
    for it = 2 : nt
        OneStepForward!(spt2, spt1, fidMtx, tmp, tmp1)
        if spt2.it <= ns
           AddSourceBorn!(spt2, path, dm, du)
        end
        CopySnapShot!(spt1, spt2)
        Spt2ShotV!(shotv, spt1)
    end
    return shotv
end

# write wavefield to hard drive, do born approximation
function MultiStepForward(path_wfd, path::String, dm::Array{Float64,2}, du::Array{Float64,2}, fidMtx::FidMtx; tmax=0.5)
    (nz, nx, ext, iflag, dt, ns) = InfoPv(path);
    nt  = round(Int64, tmax/dt)+1
    spt1= InitSnapShot(nz, nx, ext, iflag, dt, 1)
    spt2= InitSnapShot(nz, nx, ext, iflag, dt, 2)
    AddSourceBorn!(spt1, path, dm, du)
    fid = WriteWfd(path_wfd, spt1)
    tmp = zeros(length(spt1.Vxx))
    tmp1= zeros(tmp)
    for it = 2 : nt
        OneStepForward!(spt2, spt1, fidMtx, tmp, tmp1)
        if spt2.it <= ns
           AddSourceBorn!(spt2, path, dm, du)
        end
        CopySnapShot!(spt1, spt2)
        WriteWfd(fid, spt1)
    end
    close(fid)
    return nothing
end

# output shotv(Vx, Vz), input snapshot
function MultiStepForward(irz::Array{Int64,1}, irx::Array{Int64,1}, path::String, fidMtx::FidMtx; tmax=0.5)
    (nz, nx, ext, iflag, dt, ns) = InfoSnapShots(path);
    nt  = round(Int64, tmax/dt)+1
    spt1= ReadSnapShot(path, 1)
    spt2= InitSnapShot(nz, nx, ext, iflag, dt, 2)
    shotv= InitShotV(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
    Spt2ShotV!(shotv, spt1)
    tmp = zeros(length(spt1.Vxx))
    tmp1= zeros(tmp)
    for it = 2 : nt
        OneStepForward!(spt2, spt1, fidMtx, tmp, tmp1)
        if spt2.it <= ns
           spt = ReadSnapShot(path, it)
           Add2SnapShots!(spt2, spt)
        end
        CopySnapShot!(spt1, spt2)
        Spt2ShotV!(shotv, spt1)
    end
    return shotv
end

# return shot(vx, vz, Txx, Tzz, Txz), inject multiple souces
# function MultiStepForward(irz::Array{Int64,1}, irx::Array{Int64,1}, srcs::Array{Source,1}, fidMtx::FidMtx; tmax=0.5)
#     nz  = srcs[1].nz ; nx = srcs[1].nx    ;
#     ext = srcs[1].ext; iflag=srcs[1].iflag; dt = srcs[1].dt;
#     (stl, stu) = SourcesTimeRange(srcs)
#     nt  = round(Int64, tmax/dt)+1
#     spt1 = InitSnapShot(nz, nx, ext, iflag, dt, 1)
#     spt2 = InitSnapShot(nz, nx, ext, iflag, dt, 2)
#     AddMultiSources!(spt1, srcs)
#     shot = InitShot(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
#     Spt2Shot!(shot, spt1)
    # tmp = zeros(length(spt1.Vxx))
    # tmp1= zeros(tmp)
#     for it = 2 : nt
#         OneStepForward!(spt2, spt1, fidMtx, tmp, tmp1)
#         if stl <= (it-1)*dt <= stu
#            AddMultiSources!(spt2, srcs)
#         end
#         CopySnapShot!(spt1, spt2)
#         Spt2Shot!(shot, spt1)
#     end
#     return shot
# end
