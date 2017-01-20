# ========================Adjoint=========================================
# output wave field, input shotv
function MultiStepAdjoint(path::String, shotv::ShotV, fidMtx::FidMtx)
    nz=shotv.nz; nx=shotv.nx; ext=shotv.ext; iflag=shotv.iflag
    ot = shotv.ot; dt = shotv.dt; nt = shotv.nt
    spt1 = InitSnapShot(nz, nx, ext, iflag, dt, nt  )
    spt2 = InitSnapShot(nz, nx, ext, iflag, dt, nt-1)
    AddShotV2Spt!(spt1, shotv)
    fid = WriteWfd(path, spt1)
    for it = nt-1: -1: 1
        OneStepAdjoint!(spt2, spt1, fidMtx)
        AddShotV2Spt!(spt2, shotv)
        CopySnapShot!(spt1, spt2)
        WriteWfd(fid, spt1)
    end
    close(fid)
    ReverseOrderWfd(path)
    return nothing
end

function MultiStepAdjoint(shotv::ShotV, fidMtx::FidMtx, path::String)
    (nz, nx, ext, iflag, dt, ns) = InfoPv(path);
    dm = zeros(nz, nx); du = zeros(nz, nx);
    ot = shotv.ot; dt = shotv.dt; nt = shotv.nt
    spt1 = InitSnapShot(nz, nx, ext, iflag, dt, nt  )
    spt2 = InitSnapShot(nz, nx, ext, iflag, dt, nt-1)
    AddShotV2Spt!(spt1, shotv)
    if nt <= ns
       (dm, du) = imaging(dm, du, spt1, path)
    end
    for it = nt-1: -1: 1
        OneStepAdjoint!(spt2, spt1, fidMtx)
        AddShotV2Spt!(spt2, shotv)
        CopySnapShot!(spt1, spt2)
        if it <= ns
           (dm, du) = imaging(dm, du, spt1, path)
        end
    end
    return dm, du
end

# output SnapShots, input shotv
# function MultiStepAdjoint_spt(path::String, shotv::ShotV, fidMtx::FidMtx)
#     nz=shotv.nz; nx=shotv.nx; ext=shotv.ext; iflag=shotv.iflag
#     ot = shotv.ot; dt = shotv.dt; nt = shotv.nt
#     spt1 = InitSnapShot(nz, nx, ext, iflag, dt, nt  )
#     spt2 = InitSnapShot(nz, nx, ext, iflag, dt, nt-1)
#     AddShotV2Spt!(spt1, shotv)
#     fid = WriteSnapShot(path, spt1)
#     for it = nt-1: -1: 1
#         OneStepAdjoint!(spt2, spt1, fidMtx)
#         AddShotV2Spt!(spt2, shotv)
#         CopySnapShot!(spt1, spt2)
#         WriteSnapShot(fid, spt1)
#     end
#     close(fid)
#     ReverseOrderSnapShots(path)
#     return nothing
# end
