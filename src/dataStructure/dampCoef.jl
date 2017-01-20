immutable dampCoef
    dVxx  :: SparseMatrixCSC{Float64, Int64}
    dVxz  :: SparseMatrixCSC{Float64, Int64}
    dVzx  :: SparseMatrixCSC{Float64, Int64}
    dVzz  :: SparseMatrixCSC{Float64, Int64}
    dTxxx :: SparseMatrixCSC{Float64, Int64}
    dTxxz :: SparseMatrixCSC{Float64, Int64}
    dTzzx :: SparseMatrixCSC{Float64, Int64}
    dTzzz :: SparseMatrixCSC{Float64, Int64}
    dTxzx :: SparseMatrixCSC{Float64, Int64}
    dTxzz :: SparseMatrixCSC{Float64, Int64}
end

function InitDampCoef(nz::Int64, nx::Int64, ext::Int64, iflag::Int64)
  if iflag == 1
     Nz = nz + 2*ext
  elseif iflag == 2
     Nz = nz +   ext
  end
  Nx = nx + 2*ext
  dpCoef = dampCoef(spzeros(Nz,Nx), spzeros(Nz,Nx), spzeros(Nz,Nx), spzeros(Nz,Nx), spzeros(Nz,Nx),
                    spzeros(Nz,Nx), spzeros(Nz,Nx), spzeros(Nz,Nx), spzeros(Nz,Nx), spzeros(Nz,Nx))
  return dpCoef
end

function DampBound(nz::Int64, nx::Int64, maxv::Float64, ext::Int64, rho::Array{Float64,2}, iflag::Int64, dz::Float64)
    if iflag == 1
       Nz = nz + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    if ext == 5
       R = 0.01
    elseif ext == 10
       R = 0.001
    elseif ext >= 20
       R = 0.0001
    end
    ll = 2
    dVxx  = spzeros(Nz, Nx)
    dVxz  = spzeros(Nz, Nx)
    dVzx  = spzeros(Nz, Nx)
    dVzz  = spzeros(Nz, Nx)
    dTxxx = spzeros(Nz, Nx)
    dTxxz = spzeros(Nz, Nx)
    dTzzx = spzeros(Nz, Nx)
    dTzzz = spzeros(Nz, Nx)
    dTxzx = spzeros(Nz, Nx)
    dTxzz = spzeros(Nz, Nx)
    for i = 1 : ext
        tmph = -1.5 * maxv/(ext*dz) * log(R) * (((ext+0.5-i)/ext)^ll)
        tmpo = -1.5 * maxv/(ext*dz) * log(R) * (((ext+1.0-i)/ext)^ll)
        # upper
        if iflag == 1
           dVxz[i,:]  = tmpo * rho[i,:]
           dVzz[i,:]  = tmph * rho[i,:]
           dTxxz[i,:] = tmpo * ones(Nx)
           dTzzz[i,:] = tmpo * ones(Nx)
           dTxzz[i,:] = tmph * ones(Nx)
        end
        # left
        dVxx[:,i]  = tmph * rho[:,i]
        dVzx[:,i]  = tmpo * rho[:,i]
        dTxxx[:,i] = tmpo * ones(Nz)
        dTzzx[:,i] = tmpo * ones(Nz)
        dTxzx[:,i] = tmph * ones(Nz)
    end
    # lower side
    if iflag == 1
       for i = nz+ext+1 : Nz
           tmph = -1.5 * maxv/(ext*dz) * log(R) * (((i+0.5-nz-ext)/ext)^ll)
           tmpo = -1.5 * maxv/(ext*dz) * log(R) * (((i    -nz-ext)/ext)^ll)
           dVxz[i,:]  = tmpo * rho[i,:]
           dVzz[i,:]  = tmph * rho[i,:]
           dTxxz[i,:] = tmpo * ones(Nx)
           dTzzz[i,:] = tmpo * ones(Nx)
           dTxzz[i,:] = tmph * ones(Nx)
       end
    elseif iflag == 2
        for i = nz+1 : Nz
            tmph = -1.5 * maxv/(ext*dz) * log(R) * (((i+0.5-nz)/ext)^ll)
            tmpo = -1.5 * maxv/(ext*dz) * log(R) * (((i    -nz)/ext)^ll)
            dVxz[i,:]  = tmpo * rho[i,:]
            dVzz[i,:]  = tmph * rho[i,:]
            dTxxz[i,:] = tmpo * ones(Nx)
            dTzzz[i,:] = tmpo * ones(Nx)
            dTxzz[i,:] = tmph * ones(Nx)
        end
    end
    # right side
    for i = nx+ext+1 : Nx
        tmph = -1.5 * maxv/(ext*dz) * log(R) * (((i+0.5-nx-ext)/ext)^ll)
        tmpo = -1.5 * maxv/(ext*dz) * log(R) * (((i    -nx-ext)/ext)^ll)
        dVxx[:,i]  = tmph * rho[:,i]
        dVzx[:,i]  = tmpo * rho[:,i]
        dTxxx[:,i] = tmpo * ones(Nz)
        dTzzx[:,i] = tmpo * ones(Nz)
        dTxzx[:,i] = tmph * ones(Nz)
    end
    dpCoef = dampCoef(dVxx, dVxz, dVzx, dVzz, dTxxx, dTxxz, dTzzx, dTzzz, dTxzx, dTxzz)
    return dpCoef
end

function modExpand(par::Array{Float64,2}, ext::Int64, iflag::Int64)
    temp1 = repmat(par[:,  1], 1, ext)
    temp2 = repmat(par[:,end], 1, ext)
    par   = hcat(temp1, par, temp2)
    if iflag == 1
       temp1 = repmat((par[1,  :])', ext, 1)
       temp2 = repmat((par[end,:])', ext, 1)
       par   = vcat(temp1, par, temp2)
    elseif iflag == 2
       temp2 = repmat((par[end,:])', ext, 1)
       par   = vcat(par, temp2)
    end
    return par
end

function modelPadding(lambda::Array{Float64,2}, mu::Array{Float64,2}, rho::Array{Float64,2}, ext::Int64, iflag::Int64)
    lambda = modExpand(lambda, ext, iflag)
    mu     = modExpand(mu    , ext, iflag)
    rho    = modExpand(rho   , ext, iflag)
    return lambda, mu, rho
end

function modSmooth(par::Array{Float64,2}, L::Int64)
    n = 2*L + 1
    par1 = modExpand(par, L, 1)
    s = convert(Array{typeof(par1[1]), 1}, hanning(n))
    s = s / sum(s)
    par1 = conv2(s,s,par1)
    par1 = par1[2*L+1:end-2*L, 2*L+1:end-2*L]
    return par1
end

# function DampBound(nz::Int64, nx::Int64, maxv::Float64, ext::Int64, rho::Array{Float64,2}, iflag::Int64, dz::Float64)
#     if iflag == 1
#        Nz = nz + 2*ext
#     elseif iflag == 2
#        Nz = nz +   ext
#     end
#     Nx = nx + 2*ext
#     if ext == 5
#        R = 0.01
#     elseif ext == 10
#        R = 0.001
#     elseif ext >= 20
#        R = 0.0001
#     end
#     m =0.25; n = 0.75; ll = 4;
#     dVxx  = spzeros(Nz, Nx)
#     dVxz  = spzeros(Nz, Nx)
#     dVzx  = spzeros(Nz, Nx)
#     dVzz  = spzeros(Nz, Nx)
#     dTxxx = spzeros(Nz, Nx)
#     dTxxz = spzeros(Nz, Nx)
#     dTzzx = spzeros(Nz, Nx)
#     dTzzz = spzeros(Nz, Nx)
#     dTxzx = spzeros(Nz, Nx)
#     dTxzz = spzeros(Nz, Nx)
#     #  upper and left sides
#     for i = 1 : ext
#         tmph = -maxv/ext * log(R) * (m*(ext+0.5-i)/ext + n*((ext+0.5-i)/ext)^ll)
#         tmpo = -maxv/ext * log(R) * (m*(ext+1.0-i)/ext + n*((ext+1.0-i)/ext)^ll)
#         if iflag == 1
#            dVxz[i,:]  = tmpo * rho[i,:]
#            dVzz[i,:]  = tmph * rho[i,:]
#            dTxxz[i,:] = tmpo * ones(Nx)
#            dTzzz[i,:] = tmpo * ones(Nx)
#            dTxzz[i,:] = tmph * ones(Nx)
#         end
#         dVxx[:,i]  = tmph * rho[:,i]
#         dVzx[:,i]  = tmpo * rho[:,i]
#         dTxxx[:,i] = tmpo * ones(Nz)
#         dTzzx[:,i] = tmpo * ones(Nz)
#         dTxzx[:,i] = tmph * ones(Nz)
#     end
#     # lower side
#     if iflag == 1
#        for i = nz+ext+1 : Nz
#            tmph = -maxv/ext * log(R) * (m*(i+0.5-nz-ext)/ext + n*((i+0.5-nz-ext)/ext)^ll)
#            tmpo = -maxv/ext * log(R) * (m*(i    -nz-ext)/ext + n*((i    -nz-ext)/ext)^ll)
#            dVxz[i,:]  = tmpo * rho[i,:]
#            dVzz[i,:]  = tmph * rho[i,:]
#            dTxxz[i,:] = tmpo * ones(Nx)
#            dTzzz[i,:] = tmpo * ones(Nx)
#            dTxzz[i,:] = tmph * ones(Nx)
#        end
#     elseif iflag == 2
#        for i = nz+1 : Nz
#            tmph = -maxv/ext * log(R) * (m*(i+0.5-nz-ext)/ext + n*((i+0.5-nz-ext)/ext)^ll)
#            tmpo = -maxv/ext * log(R) * (m*(i    -nz-ext)/ext + n*((i    -nz-ext)/ext)^ll)
#            dVxz[i,:]  = tmpo * rho[i,:]
#            dVzz[i,:]  = tmph * rho[i,:]
#            dTxxz[i,:] = tmpo * ones(Nx)
#            dTzzz[i,:] = tmpo * ones(Nx)
#            dTxzz[i,:] = tmph * ones(Nx)
#        end
#     end
#     # right side
#     for i = nx+ext+1 : Nx
#         tmph = -maxv/ext * log(R) * (m*(i+0.5-nx-ext)/ext + n*((i+0.5-nx-ext)/ext)^ll)
#         tmpo = -maxv/ext * log(R) * (m*(i    -nx-ext)/ext + n*((i    -nx-ext)/ext)^ll)
#         dVxx[:,i]  = tmph * rho[:,i]
#         dVzx[:,i]  = tmpo * rho[:,i]
#         dTxxx[:,i] = tmpo * ones(Nz)
#         dTzzx[:,i] = tmpo * ones(Nz)
#         dTxzx[:,i] = tmph * ones(Nz)
#     end
#     dpCoef = dampCoef(dVxx, dVxz, dVzx, dVzz, dTxxx, dTxxz, dTzzx, dTzzz, dTxzx, dTxzz)
#     return dpCoef
# end
