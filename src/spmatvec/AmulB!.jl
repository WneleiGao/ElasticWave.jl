function AmulB!(y::Array{Float64,1}, A::SparseMatrixCSC{Float64,Int64}, x::Array{Float64,1})
    (m,n) = size(A)
    if length(x) != n || length(y) != m
       throw(DimensionMismatch("length(x) != n || length(y) != m"))
    end
    p = ccall((:a_mul_b_rr_, spmatveclib), Int64, (Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Float64}),
                                                   &m,         &n,         A.nzval,      A.rowval,   A.colptr,   x,            y             )
    return nothing
end

function AmulB!(y::Array{Complex128,1}, A::SparseMatrixCSC{Float64,Int64}, x::Array{Complex128,1})
    (m,n) = size(A)
    if length(x) != n || length(y) != m
       throw(DimensionMismatch("length(x) != n || length(y) != m"))
    end
    p = ccall((:a_mul_b_rc_, spmatveclib), Int64, (Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Int64}, Ptr{Int64}, Ptr{Complex128}, Ptr{Complex128}),
                                                   &m,         &n,         A.nzval,      A.rowval,   A.colptr,   x,               y                )
    return nothing
end

function AmulB!(y::Array{Complex128,1}, A::SparseMatrixCSC{Complex128,Int64}, x::Array{Complex128,1})
    (m,n) = size(A)
    if length(x) != n || length(y) != m
       throw(DimensionMismatch("length(x) != n || length(y) != m"))
    end
    p = ccall((:a_mul_b_cc_, spmatveclib), Int64, (Ptr{Int64}, Ptr{Int64}, Ptr{Complex128}, Ptr{Int64}, Ptr{Int64}, Ptr{Complex128}, Ptr{Complex128}),
                                                      &m,         &n,         A.nzval,         A.rowval,   A.colptr,   x,               y                )
    return nothing
end
