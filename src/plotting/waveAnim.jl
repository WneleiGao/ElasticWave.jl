function waveAnim(path::String, path1::String; itype="wfd", cpt="Vz",cmap="PuOr",
                  pclip=0.98, vlim=1.0, aspect="auto", interpolation="hanning", ox=0.0,
                  oy=0.0, dx=10.0, dy=10.0, dt1=0.002, wbox=8, hbox=8, spd=40)
    fid = open(path1, "r")
    nz = read(fid, Int32); nx = read(fid, Int32); nz = Int64(nz); nx = Int64(nx);
    ext= read(fid, Int32); iflag = read(fid, Int32)
    dt = read(fid, Float32);
    inc = round(Int64, dt1/dt)
    fig = plt.figure(figsize=(wbox, hbox))
    ims = PyCall.PyObject[]
    if itype == "wfd"
       sptSize = sizeof(Float32)*nz*nx*5
       nt = floor(Int64, (filesize(fid)-sizeof(Float32)*5)/sptSize)
       for it = 1 : inc : nt
           wfd = ReadWfd(path1, it)
           if cpt == "Vx"
              d = reshape(wfd.Vx, nz, nx)
           elseif cpt == "Vz"
              d = reshape(wfd.Vz, nz, nx)
           elseif cpt == "Txx"
              d = reshape(wfd.Txx, nz, nx)
           elseif cpt == "Tzz"
              d = reshape(wfd.Tzz, nz, nx)
           elseif cpt == "Txz"
              d = reshape(wfd.Txz, nz, nx)
           end
           if pclip < 1.0
              vlim = quantile(vec(abs(d)), pclip)
           end
           im = plt.imshow(d, cmap=cmap, vmin=-vlim, vmax=vlim, extent=[ox, ox+(size(d,2)-1)*dx, oy+(size(d,1)-1)*dy,oy], aspect=aspect, interpolation=interpolation)
           push!(ims, PyCall.PyObject[im])
       end
    elseif itype == "spt"
       if iflag == 1
          Nz = nz + 2*ext
          upper = ext
       elseif iflag == 2
          Nz = nz +   ext
          upper = Int32(0)
       end
       Nx  = nx + 2*ext
       sptSize = sizeof(Float32)*Nz*Nx*10
       nt = floor(Int64, (filesize(fid)-sizeof(Float32)*5)/sptSize)
       for it = 1 : inc : nt
           spt = ReadSnapShot(path1, it)
           if cpt == "Vx"
              tmp = reshape(spt.Vxx+spt.Vxz, Nz, Nx)
              d = tmp[upper+1:upper+nz, ext+1:nx+ext]
           elseif cpt == "Vz"
              tmp = reshape(spt.Vzx+spt.Vzz, Nz, Nx)
              d = tmp[upper+1:upper+nz, ext+1:nx+ext]
           elseif cpt == "Txx"
              tmp = reshape(spt.Txxx+spt.Txxz, Nz, Nx)
              d = tmp[upper+1:upper+nz, ext+1:nx+ext]
           elseif cpt == "Tzz"
              tmp = reshape(spt.Tzzx+spt.Tzzz, Nz, Nx)
              d = tmp[upper+1:upper+nz, ext+1:nx+ext]
           elseif cpt == "Txz"
              tmp = reshape(spt.Txzx+spt.Txzz, Nz, Nx)
              d = tmp[upper+1:upper+nz, ext+1:nx+ext]
           end
           if pclip < 1.0
              vlim = quantile(vec(abs(d)), pclip)
           end
           im = plt.imshow(d, cmap=cmap, vmin=-vlim, vmax=vlim, extent=[ox, ox+(size(d,2)-1)*dx, oy+(size(d,1)-1)*dy,oy], aspect=aspect, interpolation=interpolation)
           push!(ims, PyCall.PyObject[im])
       end
    end
    ani = anim.ArtistAnimation(fig, ims, interval=spd, blit=true)
    path = join([path "_" cpt ".mp4"])
    ani[:save](path)
end
