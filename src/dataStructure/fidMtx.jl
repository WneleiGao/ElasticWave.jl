immutable FidMtx
  MVxxBVxx   :: Array{Float64,1}
  MVxxBTxx   :: SparseMatrixCSC{Float64,Int64}
  MVxzBVxz   :: Array{Float64,1}
  MVxzBTxz   :: SparseMatrixCSC{Float64,Int64}
  MVzxBVzx   :: Array{Float64,1}
  MVzxBTxz   :: SparseMatrixCSC{Float64,Int64}
  MVzzBVzz   :: Array{Float64,1}
  MVzzBTzz   :: SparseMatrixCSC{Float64,Int64}
  MTxxxBVx   :: SparseMatrixCSC{Float64,Int64}
  MTxxxBTxxx :: Array{Float64,1}
  MTxxzBVz   :: SparseMatrixCSC{Float64,Int64}
  MTxxzBTxxz :: Array{Float64,1}
  MTzzxBVx   :: SparseMatrixCSC{Float64,Int64}
  MTzzxBTzzx :: Array{Float64,1}
  MTzzzBVz   :: SparseMatrixCSC{Float64,Int64}
  MTzzzBTzzz :: Array{Float64,1}
  MTxzxBVz   :: SparseMatrixCSC{Float64,Int64}
  MTxzxBTxzx :: Array{Float64,1}
  MTxzzBVx   :: SparseMatrixCSC{Float64,Int64}
  MTxzzBTxzz :: Array{Float64,1}
end

function vel2lm(vp::Array{Float64,2}, vs::Array{Float64,2}, rho::Array{Float64,2})
    (m, n) = size(vp)
    lambda = zeros(typeof(vp[1]), m, n)
    mu     = zeros(typeof(vp[1]), m, n)
    for j = 1: n
        for i = 1: m
            mu[i,j] = rho[i,j] * vs[i,j]^2
            lambda[i,j] = rho[i,j] * vp[i,j]^2 - 2*mu[i,j]
        end
    end
    return lambda, mu
end

function lm2vel(lambda::Array{Float64,2}, mu::Array{Float64,2}, rho::Array{Float64,2})
    (m, n) = size(lambda)
    vp = Array(typeof(lambda[1]), m, n)
    vs = Array(typeof(lambda[1]), m, n)
    for j = 1: n
        for i = 1: m
            vp[i,j] = sqrt((lambda[i,j]+2*mu[i,j])/rho[i,j])
            vs[i,j] = sqrt(mu[i,j]/rho[i,j])
        end
    end
    return vp, vs
end

function DispStable!(vmax::Float64, vmin::Float64, f0::Float64, dz::Float64, dx::Float64, dt::Float64)
    h = vmin / 5 / f0
    d = minimum([dz dx])
    tt = 6*d / (7*sqrt(3)*vmax)
    if dx >= h || dz >= h || dt >= tt
       println("maximum acceptable grid size: $h")
       println("maximum acceptable time step: $tt")
       error("unstable or frequency dispersion")
    end
    return nothing
end

function MVxxBVxx(dVxx::SparseMatrixCSC{Float64,Int64}, rho::Array{Float64,2}, dt::Float64)
    (m,n) = size(dVxx)
    Mdiag = zeros(m*n)
    for ix = 1 : n-1
        for iz = 1 : m
            a3 = (rho[iz,ix+1]+rho[iz,ix])/2/dt + dVxx[iz,ix]/2
            tmp= (rho[iz,ix+1]+rho[iz,ix])/2/(dt*a3) - dVxx[iz,ix]/(2*a3)
            Mdiag[(ix-1)*m+iz] = tmp
        end
    end
    for iz = 1 : m
        a3 = rho[iz,n]/dt + dVxx[iz,n]/2
        tmp= rho[iz,n]/(dt*a3) - dVxx[iz,n]/(2*a3)
        Mdiag[(n-1)*m+iz] = tmp
    end
    M = Mdiag
    return M
end

function MVxxBTxx(dVxx::SparseMatrixCSC{Float64,Int64}, rho::Array{Float64,2}, dt::Float64, dx::Float64)
    a1 =  9/8
    a2 = -1/24
    (m, n) = size(dVxx)
    A3 = zeros(m*n)
    for ix = 1:n-1
        for iz = 1:m
            A3[(ix-1)*m+iz] = 1 / ((rho[iz,ix+1]+rho[iz,ix])/2/dt + dVxx[iz,ix]/2)
        end
    end
    for iz = 1 : m
        A3[(n-1)*m+iz] =  1 / (rho[iz,n]/dt + dVxx[iz,n]/2)
    end
    A = spzeros(n, n)
    A[1, 1] = -a1/dx; A[1, 2] = a1/dx; A[1, 3] = a2/dx;
    for ix = 2 : n-2
        A[ix, ix-1] = -a2/dx
        A[ix, ix  ] = -a1/dx
        A[ix, ix+1] =  a1/dx
        A[ix, ix+2] =  a2/dx
    end
    A[n-1, n-2] = -a2/dx;
    A[n-1, n-1] = -a1/dx; A[n-1, n  ] = a1/dx;
    A[n  , n-1] = -a2/dx; A[n  , n  ] = -a1/dx;
    Im = speye(m)
    M = kron(A, Im)
    M = spdiagm(A3) * M
    return M
end

function MVxzBVxz(dVxz::SparseMatrixCSC{Float64,Int64}, rho::Array{Float64,2}, dt::Float64)
    (m,n) = size(dVxz)
    Mdiag = zeros(m*n)
    for ix = 1 : n-1
        for iz = 1 : m
            a3 = (rho[iz,ix+1]+rho[iz,ix])/2/dt + dVxz[iz,ix]/2
            tmp= (rho[iz,ix+1]+rho[iz,ix])/2/(dt*a3) - dVxz[iz,ix]/(2*a3)
            Mdiag[(ix-1)*m+iz] = tmp
        end
    end
    for iz = 1 : m
        a3 = rho[iz,n]/dt + dVxz[iz,n]/2
        tmp= rho[iz,n]/(dt*a3) - dVxz[iz,n]/(2*a3)
        Mdiag[(n-1)*m+iz] = tmp
    end
    M = Mdiag
    return M
end

function MVxzBTxz(dVxz::SparseMatrixCSC{Float64, Int64}, rho::Array{Float64,2}, dt::Float64, dz::Float64)
    a1 =  9/8
    a2 = -1/24
    (m,n) = size(dVxz)
    A3 = zeros(m*n)
    for ix = 1:n-1
        for iz = 1:m
            A3[(ix-1)*m+iz] = 1 / ((rho[iz,ix+1]+rho[iz,ix])/2/dt + dVxz[iz,ix]/2)
        end
    end
    for iz = 1 : m
        A3[(n-1)*m+iz] =  1 / (rho[iz,n]/dt + dVxz[iz,n]/2)
    end
    A = spzeros(m, m)
    A[1, 1] = a1/dz; A[1, 2] = a2/dz;
    A[2, 1] =-a1/dz; A[2, 2] = a1/dz; A[2, 3] = a2/dz;
    for iz = 3 : m-1
        A[iz, iz-2] = -a2/dz
        A[iz, iz-1] = -a1/dz
        A[iz, iz  ] =  a1/dz
        A[iz, iz+1] =  a2/dz
    end
    A[m, m-2] = -a2/dz;
    A[m, m-1] = -a1/dz; A[m, m] = a1/dz;
    In = speye(n)
    M = kron(In, A)
    M = spdiagm(A3) * M
    return M
end

function MVzxBVzx(dVzx::SparseMatrixCSC{Float64,Int64}, rho::Array{Float64,2}, dt::Float64)
    (m,n) = size(dVzx)
    Mdiag = zeros(m*n)
    for ix = 1 : n
        for iz = 1 : m-1
            a3 = (rho[iz+1,ix]+rho[iz,ix])/2/dt + dVzx[iz,ix]/2
            tmp= (rho[iz+1,ix]+rho[iz,ix])/2/(dt*a3) - dVzx[iz,ix]/(2*a3)
            Mdiag[(ix-1)*m+iz] = tmp
        end
    end
    for ix = 1 : n
        a3 = rho[m,ix]/dt + dVzx[m,ix]/2
        tmp= rho[m,ix]/(dt*a3) - dVzx[m,ix]/(2*a3)
        Mdiag[(ix-1)*m+m] = tmp
    end
    M = Mdiag
    return M
end

function MVzxBTxz(dVzx::SparseMatrixCSC{Float64, Int64}, rho::Array{Float64,2}, dt::Float64, dx::Float64)
    a1 = 9/8
    a2 = -1/24
    (m,n) = size(dVzx)
    A3 = zeros(m*n)
    for ix = 1:n
        for iz = 1:m-1
            A3[(ix-1)*m+iz] = 1 / ((rho[iz+1,ix]+rho[iz,ix])/2/dt+dVzx[iz,ix]/2)
        end
    end
    for ix = 1 : n
        A3[(ix-1)*m+m] =  1 / (rho[m,ix]/dt + dVzx[m,ix]/2)
    end
    A = spzeros(n, n)
    A[1, 1] = a1/dx; A[1, 2] = a2/dx;
    A[2, 1] =-a1/dx; A[2, 2] = a1/dx;
    A[2, 3] = a2/dx;
    for ix = 3:n-1
        A[ix, ix-2] = -a2/dx
        A[ix, ix-1] = -a1/dx
        A[ix, ix  ] =  a1/dx
        A[ix, ix+1] =  a2/dx
    end
    A[n, n-2] =-a2/dx;
    A[n, n-1] =-a1/dx; A[n, n] = a1/dx;
    Im = speye(m)
    M = kron(A, Im)
    M = spdiagm(A3) * M
    return M
end

function MVzzBVzz(dVzz::SparseMatrixCSC{Float64,Int64}, rho::Array{Float64,2}, dt::Float64)
    (m,n) = size(dVzz)
    Mdiag = zeros(m*n)
    for ix = 1 : n
        for iz = 1 : m-1
            a3 = (rho[iz+1,ix]+rho[iz,ix])/2/dt + dVzz[iz,ix]/2
            tmp= (rho[iz+1,ix]+rho[iz,ix])/2/(dt*a3) - dVzz[iz,ix]/(2*a3)
            Mdiag[(ix-1)*m+iz] = tmp
        end
    end
    for ix = 1 : n
        a3 = rho[m,ix]/dt + dVzz[m,ix]/2
        tmp= rho[m,ix]/(dt*a3) - dVzz[m,ix]/(2*a3)
        Mdiag[(ix-1)*m+m] = tmp
    end
    M = Mdiag
    return M
end

function MVzzBTzz(dVzz::SparseMatrixCSC{Float64,Int64}, rho::Array{Float64,2}, dt::Float64, dz::Float64)
    a1 = 9/8
    a2 = -1/24
    (m,n) = size(dVzz)
    A3 = zeros(m*n)
    for ix = 1:n
        for iz = 1:m-1
            A3[(ix-1)*m+iz] = 1 / ((rho[iz+1,ix]+rho[iz,ix])/2/dt + dVzz[iz,ix]/2)
        end
    end
    for ix = 1 : n
        A3[(ix-1)*m+m] =  1 / (rho[m,ix]/dt + dVzz[m,ix]/2)
    end
    A = spzeros(m, m)
    A[1, 1] =-a1/dz; A[1, 2] = a1/dz; A[1, 3] = a2/dz;
    for iz = 2 : m-2
        A[iz, iz-1] = -a2/dz
        A[iz, iz  ] = -a1/dz
        A[iz, iz+1] =  a1/dz
        A[iz, iz+2] =  a2/dz
    end
    A[m-1, m-2] = -a2/dz;
    A[m-1, m-1] = -a1/dz; A[m-1, m] = a1/dz;
    A[m  , m-1] = -a2/dz; A[m  , m] = -a1/dz;
    In = speye(n)
    M = kron(In, A)
    M = spdiagm(A3) * M
    return M
end

function MTxxxBTxxx(dTxxx::SparseMatrixCSC{Float64,Int64}, dt::Float64)
    (m,n) = size(dTxxx)
    Mdiag = zeros(m*n)
    for ix = 1 : n
        for iz = 1 : m
            a3 = 1/dt + dTxxx[iz,ix]/2
            tmp= 1/(dt*a3) - dTxxx[iz,ix]/(2*a3)
            Mdiag[(ix-1)*m+iz] = tmp
        end
    end
    M = Mdiag
    return M
end

function MTxxxBVx(dTxxx::SparseMatrixCSC{Float64,Int64}, lambda::Array{Float64,2}, mu::Array{Float64,2}, dt::Float64, dx::Float64)
    a1 = 9/8
    a2 = -1/24
    (m,n) = size(dTxxx)
    A3 = zeros(m*n)
    for ix = 1 : n
        for iz = 1 : m
            A3[(ix-1)*m+iz] = (lambda[iz,ix]+2*mu[iz,ix]) / (1/dt+dTxxx[iz,ix]/2)
        end
    end
    A = spzeros(n, n)
    A[1, 1] = a1/dx; A[1, 2] = a2/dx;
    A[2, 1] =-a1/dx; A[2, 2] = a1/dx; A[2, 3] = a2/dx;
    for ix = 3:n-1
        A[ix, ix-2] = -a2/dx
        A[ix, ix-1] = -a1/dx
        A[ix, ix  ] =  a1/dx
        A[ix, ix+1] =  a2/dx
    end
    A[n, n-2] = -a2/dx;
    A[n, n-1] = -a1/dx; A[n, n] = a1/dx;
    Im = speye(m)
    M = kron(A, Im)
    M = spdiagm(A3) * M
    return M
end

function MTxxzBTxxz(dTxxz::SparseMatrixCSC{Float64,Int64}, dt::Float64)
    (m,n) = size(dTxxz)
    Mdiag = zeros(m*n)
    for ix = 1:n
        for iz = 1:m
            a3 = 1/dt + dTxxz[iz,ix]/2
            tmp= 1/(dt*a3) - dTxxz[iz,ix]/(2*a3)
            Mdiag[(ix-1)*m+iz] = tmp
        end
    end
    M = Mdiag
    return M
end

function MTxxzBVz(dTxxz::SparseMatrixCSC{Float64,Int64}, lambda::Array{Float64,2}, dt::Float64, dz::Float64)
    a1 = 9/8
    a2 = -1/24
    (m,n) = size(dTxxz)
    A3 = zeros(m*n)
    for ix = 1:n
        for iz = 1:m
            A3[(ix-1)*m+iz] = lambda[iz,ix] / (1/dt+dTxxz[iz,ix]/2)
        end
    end
    A = spzeros(m, m)
    A[1, 1] = a1/dz; A[1, 2] = a2/dz;
    A[2, 1] =-a1/dz; A[2, 2] = a1/dz; A[2, 3] = a2/dz;
    for iz = 3:m-1
        A[iz, iz-2] = -a2/dz
        A[iz, iz-1] = -a1/dz
        A[iz, iz  ] =  a1/dz
        A[iz, iz+1] =  a2/dz
    end
    A[m, m-2] = -a2/dz;
    A[m, m-1] = -a1/dz; A[m, m] = a1/dz;
    In = speye(n)
    M = kron(In, A)
    M = spdiagm(A3) * M
    return M
end

function MTzzxBTzzx(dTzzx::SparseMatrixCSC{Float64,Int64}, dt::Float64)
    (m,n) = size(dTzzx)
    Mdiag = zeros(m*n)
    for ix = 1 : n
        for iz = 1 : m
            a3 = 1/dt + dTzzx[iz,ix]/2
            tmp= 1/(dt*a3) - dTzzx[iz,ix]/(2*a3)
            Mdiag[(ix-1)*m+iz] = tmp
        end
    end
    M = Mdiag
    return M
end

function MTzzxBVx(dTzzx::SparseMatrixCSC{Float64,Int64}, lambda::Array{Float64,2}, dt::Float64, dx::Float64)
    a1 = 9/8
    a2 =-1/24
    (m, n) = size(dTzzx)
    A3 = zeros(m*n)
    for ix = 1:n
        for iz = 1:m
            A3[(ix-1)*m+iz] = lambda[iz,ix] / (1/dt+dTzzx[iz,ix]/2)
        end
    end
    A = spzeros(n, n)
    A[1, 1] = a1/dx; A[1, 2] = a2/dx;
    A[2, 1] =-a1/dx; A[2, 2] = a1/dx; A[2, 3] = a2/dx;
    for ix = 3:n-1
        A[ix, ix-2] = -a2/dx
        A[ix, ix-1] = -a1/dx
        A[ix, ix  ] =  a1/dx
        A[ix, ix+1] =  a2/dx
    end
    A[n, n-2] = -a2/dx;
    A[n, n-1] = -a1/dx; A[n, n] = a1/dx;
    Im = speye(m)
    M = kron(A, Im)
    M = spdiagm(A3) * M
    return M
end

function MTzzzBTzzz(dTzzz::SparseMatrixCSC{Float64,Int64}, dt::Float64)
    (m,n) = size(dTzzz)
    Mdiag = zeros(m*n)
    for ix = 1:n
        for iz = 1:m
            a3 = 1/dt + dTzzz[iz,ix]/2
            tmp= 1/(dt*a3) - dTzzz[iz,ix]/(2*a3)
            Mdiag[(ix-1)*m+iz] = tmp
        end
    end
    M = Mdiag
    return M
end

function MTzzzBVz(dTzzz::SparseMatrixCSC{Float64,Int64}, lambda::Array{Float64,2}, mu::Array{Float64,2}, dt::Float64, dz::Float64)
    a1 = 9/8
    a2 = -1/24
    (m,n) = size(dTzzz)
    A3 = zeros(m*n)
    for ix = 1:n
        for iz = 1:m
            A3[(ix-1)*m+iz] = (lambda[iz,ix]+2*mu[iz,ix]) / (1/dt+dTzzz[iz,ix]/2)
        end
    end
    A = spzeros(m, m)
    A[1, 1] = a1/dz; A[1, 2] = a2/dz;
    A[2, 1] =-a1/dz; A[2, 2] = a1/dz; A[2, 3] = a2/dz;
    for iz = 3 : m-1
        A[iz, iz-2] = -a2/dz
        A[iz, iz-1] = -a1/dz
        A[iz, iz  ] =  a1/dz
        A[iz, iz+1] =  a2/dz
    end
    A[m, m-2] = -a2/dz;
    A[m, m-1] = -a1/dz; A[m, m] = a1/dz;
    In = speye(n)
    M = kron(In, A)
    M = spdiagm(A3) * M
    return M
end

function MTxzxBTxzx(dTxzx::SparseMatrixCSC{Float64,Int64}, dt::Float64)
    (m,n) = size(dTxzx)
    Mdiag = zeros(m*n)
    for ix = 1 : n
        for iz = 1 : m
            a3 = 1/dt + dTxzx[iz,ix]/2
            tmp= 1/(dt*a3) - dTxzx[iz,ix]/(2*a3)
            Mdiag[(ix-1)*m+iz] = tmp
        end
    end
    M = Mdiag
    return M
end

function MTxzxBVz(dTxzx::SparseMatrixCSC{Float64,Int64}, mu::Array{Float64,2}, dt::Float64, dx::Float64)
    a1 = 9/8
    a2 =-1/24
    (m, n) = size(dTxzx)
    A3 = zeros(m*n)
    for ix = 1 : n-1
        for iz = 1 : m-1
            A3[(ix-1)*m+iz] = (mu[iz,ix]+mu[iz+1,ix]+mu[iz,ix+1]+mu[iz+1,ix+1])/4 / (1/dt+dTxzx[iz,ix]/2)
        end
    end
    for iz = 1 : m-1
        A3[(n-1)*m+iz] = (mu[iz,n]+mu[iz+1,n])/2 / (1/dt+dTxzx[iz,n]/2)
    end
    for ix = 1 : n-1
        A3[(ix-1)*m+m] = (mu[m,ix]+mu[m,ix+1])/2 / (1/dt+dTxzx[m,ix]/2)
    end
    A3[end] = mu[m,n] / (1/dt+dTxzx[m,n]/2)
    A = spzeros(n, n)
    A[1,1] =-a1/dx; A[1,2] = a1/dx; A[1,3] = a2/dx;
    for ix = 2 : n-2
        A[ix, ix-1] = -a2/dx
        A[ix, ix  ] = -a1/dx
        A[ix, ix+1] =  a1/dx
        A[ix, ix+2] =  a2/dx
    end
    A[n-1,n-2] =-a2/dx;
    A[n-1,n-1] =-a1/dx; A[n-1,n] = a1/dx;
    A[n  ,n-1] =-a2/dx; A[n  ,n] =-a1/dx;
    Im = speye(m)
    M = kron(A, Im)
    M = spdiagm(A3) * M
    return M
end

function MTxzzBTxzz(dTxzz::SparseMatrixCSC{Float64,Int64}, dt::Float64)
    (m, n) = size(dTxzz)
    Mdiag  = zeros(m*n)
    for ix = 1:n
        for iz = 1:m
            a3 = 1/dt + dTxzz[iz,ix]/2
            tmp= 1/(dt*a3) - dTxzz[iz,ix]/(2*a3)
            Mdiag[(ix-1)*m+iz] = tmp
        end
    end
    M = Mdiag
    return M
end

function MTxzzBVx(dTxzz::SparseMatrixCSC{Float64,Int64}, mu::Array{Float64,2}, dt::Float64, dz::Float64)
    a1 = 9/8
    a2 =-1/24
    (m,n) = size(dTxzz)
    A3 = zeros(m*n)
    for ix = 1 : n-1
        for iz = 1 : m-1
            A3[(ix-1)*m+iz] = (mu[iz,ix]+mu[iz+1,ix]+mu[iz,ix+1]+mu[iz+1,ix+1])/4 / (1/dt+dTxzz[iz,ix]/2)
        end
    end
    for iz = 1 : m-1
        A3[(n-1)*m+iz] = (mu[iz,n]+mu[iz+1,n])/2 / (1/dt+dTxzz[iz,n]/2)
    end
    for ix = 1 : n-1
        A3[(ix-1)*m+m] = (mu[m,ix]+mu[m,ix+1])/2 / (1/dt+dTxzz[m,ix]/2)
    end
    A3[end] = mu[m,n] / (1/dt+dTxzz[m,n]/2)
    A = spzeros(m, m)
    A[1,1] =-a1/dz; A[1,2] = a1/dz; A[1,3] = a2/dz;
    for iz = 2:m-2
        A[iz, iz-1] = -a2/dz
        A[iz, iz  ] = -a1/dz
        A[iz, iz+1] =  a1/dz
        A[iz, iz+2] =  a2/dz
    end
    A[m-1,m-2] =-a2/dz;
    A[m-1,m-1] =-a1/dz; A[m-1,m] = a1/dz;
    A[m  ,m-1] =-a2/dz; A[m  ,m] =-a1/dz;
    In = speye(n)
    M = kron(In, A)
    M = spdiagm(A3) * M
    return M
end

function CreateFidMtx(nz::Int64, nx::Int64, ext::Int64, iflag::Int64, vp::Array{Float64,2}, vs::Array{Float64,2}, rho::Array{Float64,2}, dz::Float64, dx::Float64, dt::Float64, f0::Float64)
    vmax = maximum(vp); vmin = minimum(vs);
    DispStable!(vmax, vmin, f0, dz, dx, dt)
    (lambda, mu) = vel2lm(vp, vs, rho)
    (lambda, mu, rho) = modelPadding(lambda, mu, rho, ext, iflag)
    dpCoef = DampBound(nz, nx, vmax, ext, rho, iflag, dz)
    MVxxBVxx1 = MVxxBVxx(dpCoef.dVxx, rho, dt)
    MVxxBTxx1 = MVxxBTxx(dpCoef.dVxx, rho, dt, dx)
    MVxzBVxz1 = MVxzBVxz(dpCoef.dVxz, rho, dt)
    MVxzBTxz1 = MVxzBTxz(dpCoef.dVxz, rho, dt, dz)
    MVzxBVzx1 = MVzxBVzx(dpCoef.dVzx, rho, dt)
    MVzxBTxz1 = MVzxBTxz(dpCoef.dVzx, rho, dt, dx)
    MVzzBVzz1 = MVzzBVzz(dpCoef.dVzz, rho, dt)
    MVzzBTzz1 = MVzzBTzz(dpCoef.dVzz, rho, dt, dz)
    MTxxxBVx1   =   MTxxxBVx(dpCoef.dTxxx, lambda, mu, dt, dx)
    MTxxxBTxxx1 = MTxxxBTxxx(dpCoef.dTxxx, dt)
    MTxxzBVz1   =   MTxxzBVz(dpCoef.dTxxz, lambda, dt, dz)
    MTxxzBTxxz1 = MTxxzBTxxz(dpCoef.dTxxz, dt)
    MTzzxBVx1   =   MTzzxBVx(dpCoef.dTzzx, lambda, dt, dx)
    MTzzxBTzzx1 = MTzzxBTzzx(dpCoef.dTzzx, dt)
    MTzzzBVz1   =   MTzzzBVz(dpCoef.dTzzz, lambda, mu, dt, dz)
    MTzzzBTzzz1 = MTzzzBTzzz(dpCoef.dTzzz, dt)
    MTxzxBVz1   =   MTxzxBVz(dpCoef.dTxzx, mu, dt, dx)
    MTxzxBTxzx1 = MTxzxBTxzx(dpCoef.dTxzx, dt)
    MTxzzBVx1   =   MTxzzBVx(dpCoef.dTxzz, mu, dt, dz)
    MTxzzBTxzz1 = MTxzzBTxzz(dpCoef.dTxzz, dt)
    fidMtx = FidMtx(MVxxBVxx1, MVxxBTxx1,   MVxzBVxz1, MVxzBTxz1,
                    MVzxBVzx1, MVzxBTxz1,   MVzzBVzz1, MVzzBTzz1,
                    MTxxxBVx1, MTxxxBTxxx1, MTxxzBVz1, MTxxzBTxxz1,
                    MTzzxBVx1, MTzzxBTzzx1, MTzzzBVz1, MTzzzBTzzz1,
                    MTxzxBVz1, MTxzxBTxzx1, MTxzzBVx1, MTxzzBTxzz1)
    return fidMtx
end
