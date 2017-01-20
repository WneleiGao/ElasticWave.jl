function WriteImage(path, dm, du)
    (m, n) = size(dm)
    fid = open(path, "w")
    write(fid, Int32(m), Int32(n))
    write(fid, convert(Array{Float32},vec(dm)))
    write(fid, convert(Array{Float32},vec(du)))
    close(fid)
end

function ReadImage(path)
    fid = open(path, "r")
    m = Int64(read(fid, Int32))
    n = Int64(read(fid, Int32))
    dm = read(fid, Float32, m*n)
    dm = convert(Array{Float64}, dm)
    du = read(fid, Float32, m*n)
    du = convert(Array{Float64}, du)
    dm = reshape(dm, m, n)
    du = reshape(du, m, n)
    close(fid)
    return dm, du
end

function laplaceFilter(I)
    (m, n) = size(I)
    d = zeros(typeof(I[1]), m, n)
    for ix = 2 : n-1
        for iz = 2 : m-1
            d[iz, ix] = I[iz-1,ix]+I[iz+1,ix]+I[iz,ix-1]+I[iz,ix+1]-4*I[iz,ix]
        end
    end
    return d
end

function lambdamu2vpvs(dm, du, vps, vss, rho)
    (m, n) = size(dm)
    dvp = zeros(typeof(dm[1]), m, n)
    dvp = 2* vps .* rho .* dm
    dvs = -4 * vss .* rho .* dm + 2 * vss .* rho .* du
    return dvp, dvs
end

function envelope(d)
    (nt, nx) = size(d)
    dout     = zeros(nt, nx)
    for ix =1 : nx
        temp = hilbert(d[:, ix])
        dout[:, ix] = abs(temp)
    end
    return dout
end

function gec(din::Array{Float64,2})
    (nt, nx) = size(din)
    dout     = zeros(nt, nx)
    nop      = 2*int(nt/10) + 1
    smoother = triang(nop)
    trip     = int(nop / 20)
    env = envelope(din)
    for ix =1 : nx
        temp = env[:, ix]
        temp = conv(temp, smoother)
        envsm= temp[int(round(nop/2)) : nt+int(round(nop)/2)-1]
        envsm[1:trip] = envsm[trip+1]
        envsm[nt-trip+1: end] = envsm[nt-trip]
        trin  = din[:, ix]
        trout = trin ./ envsm
        trout = trout / maximum(abs(trout))
        dout[:, ix] = trout
    end
    return dout
end
