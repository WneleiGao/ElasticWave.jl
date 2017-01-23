function AcmulB!(y::Array{Float64,1}, A::SparseMatrixCSC{Float64,Int64}, x::Array{Float64,1}; nthreads::Int64=4)

    (m,n) = size(A)
    if length(x) != m || length(y) != n
       throw(DimensionMismatch("length(x) != m || length(y) != n"))
    end
    p = ccall((:ac_mul_b_rr_, spmatveclib), Int64, (Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Float64}),
                                                    &nthreads,  &m,         &n,         A.nzval,      A.rowval,   A.colptr,   x,            y             )
    return nothing
end

function AcmulB!(y::Array{Complex128,1}, A::SparseMatrixCSC{Float64,Int64}, x::Array{Complex128,1}; nthreads::Int64=4)

    (m,n) = size(A)
    if length(x) != m || length(y) != n
       throw(DimensionMismatch("length(x) != m || length(y) != n"))
    end
    p = ccall((:ac_mul_b_rc_, spmatveclib), Int64, (Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Int64}, Ptr{Int64}, Ptr{Complex128}, Ptr{Complex128}),
                                                       &nthreads,  &m,         &n,         A.nzval,      A.rowval,   A.colptr,   x,               y                )
    return nothing
end

function AcmulB!(y::Array{Complex128,1}, A::SparseMatrixCSC{Complex128,Int64}, x::Array{Complex128,1}; nthreads::Int64=4)

    (m,n) = size(A)
    if length(x) != m || length(y) != n
       throw(DimensionMismatch("length(x) != m || length(y) != n"))
    end
    p = ccall((:ac_mul_b_cc_, spmatveclib), Int64, (Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Complex128}, Ptr{Int64}, Ptr{Int64}, Ptr{Complex128}, Ptr{Complex128}),
                                                       &nthreads,  &m,         &n,         A.nzval,         A.rowval,   A.colptr,   x,               y                )
    return nothing
end
