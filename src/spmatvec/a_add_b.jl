function ca_add_b!(c::Array{Float64,1}, a::Array{Float64,1}, b::Array{Float64,1}; nthreads::Int64=4)
    n = length(a)
    p = ccall((:ca_add_b_rr_, spmatveclib), Int64, (Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
                                                    &nthreads , &n        , c           , a           , b             )
end
