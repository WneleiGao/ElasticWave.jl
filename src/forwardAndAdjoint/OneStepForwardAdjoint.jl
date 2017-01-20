function dotproduct!(y::Array{Float64,1}, a::Array{Float64,1}, b::Array{Float64,1})
    n = length(y)
    @inbounds for i = 1 : n
       y[i] = a[i] * b[i]
    end
    return nothing
end

function mysum1!(y::Array{Float64,1}, a::Array{Float64,1}, b::Array{Float64,1})
    n = length(y)
    @inbounds for i = 1 : n
       y[i] = a[i] + b[i]
    end
    return nothing
end

function mysum2!(a::Array{Float64,1}, b::Array{Float64,1})
    n = length(a)
    @inbounds for i = 1 : n
       a[i] = a[i] + b[i]
    end
    return nothing
end

function OneStepForward!(snapShot2::SnapShot, snapShot1::SnapShot, fidMtx::FidMtx, tmp::Array{Float64,1}, tmp1::Array{Float64,1})
    snapShot2.it  = snapShot1.it + 1
    dotproduct!(snapShot2.Vxx, fidMtx.MVxxBVxx, snapShot1.Vxx); mysum1!(tmp, snapShot1.Txxx, snapShot1.Txxz);
    AmulB!(tmp1, fidMtx.MVxxBTxx, tmp);                         mysum2!(snapShot2.Vxx, tmp1);

    mysum1!(tmp, snapShot1.Txzx, snapShot1.Txzz);
    dotproduct!(snapShot2.Vxz, fidMtx.MVxzBVxz, snapShot1.Vxz); AmulB!(tmp1, fidMtx.MVxzBTxz, tmp); mysum2!(snapShot2.Vxz, tmp1);
    dotproduct!(snapShot2.Vzx, fidMtx.MVzxBVzx, snapShot1.Vzx); AmulB!(tmp1, fidMtx.MVzxBTxz, tmp); mysum2!(snapShot2.Vzx, tmp1);

    dotproduct!(snapShot2.Vzz, fidMtx.MVzzBVzz, snapShot1.Vzz); mysum1!(tmp, snapShot1.Tzzx, snapShot1.Tzzz);
    AmulB!(tmp1, fidMtx.MVzzBTzz, tmp);                         mysum2!(snapShot2.Vzz, tmp1);

    mysum1!(tmp, snapShot2.Vxx, snapShot2.Vxz);
    dotproduct!(snapShot2.Txxx, fidMtx.MTxxxBTxxx, snapShot1.Txxx); AmulB!(tmp1, fidMtx.MTxxxBVx, tmp); mysum2!(snapShot2.Txxx, tmp1);
    dotproduct!(snapShot2.Tzzx, fidMtx.MTzzxBTzzx, snapShot1.Tzzx); AmulB!(tmp1, fidMtx.MTzzxBVx, tmp); mysum2!(snapShot2.Tzzx, tmp1);
    dotproduct!(snapShot2.Txzz, fidMtx.MTxzzBTxzz, snapShot1.Txzz); AmulB!(tmp1, fidMtx.MTxzzBVx, tmp); mysum2!(snapShot2.Txzz, tmp1);

    mysum1!(tmp, snapShot2.Vzx, snapShot2.Vzz);
    dotproduct!(snapShot2.Txxz, fidMtx.MTxxzBTxxz, snapShot1.Txxz); AmulB!(tmp1, fidMtx.MTxxzBVz, tmp); mysum2!(snapShot2.Txxz, tmp1);
    dotproduct!(snapShot2.Tzzz, fidMtx.MTzzzBTzzz, snapShot1.Tzzz); AmulB!(tmp1, fidMtx.MTzzzBVz, tmp); mysum2!(snapShot2.Tzzz, tmp1);
    dotproduct!(snapShot2.Txzx, fidMtx.MTxzxBTxzx, snapShot1.Txzx); AmulB!(tmp1, fidMtx.MTxzxBVz, tmp); mysum2!(snapShot2.Txzx, tmp1);
    return nothing
end

function OneStepForward!(snapShot2::SnapShot, snapShot1::SnapShot, fidMtx::FidMtx)
    snapShot2.it  = snapShot1.it + 1
    snapShot2.Vxx[:] = fidMtx.MVxxBVxx .* snapShot1.Vxx + fidMtx.MVxxBTxx*(snapShot1.Txxx+snapShot1.Txxz)
    snapShot2.Vxz[:] = fidMtx.MVxzBVxz .* snapShot1.Vxz + fidMtx.MVxzBTxz*(snapShot1.Txzx+snapShot1.Txzz)
    snapShot2.Vzx[:] = fidMtx.MVzxBVzx .* snapShot1.Vzx + fidMtx.MVzxBTxz*(snapShot1.Txzx+snapShot1.Txzz)
    snapShot2.Vzz[:] = fidMtx.MVzzBVzz .* snapShot1.Vzz + fidMtx.MVzzBTzz*(snapShot1.Tzzx+snapShot1.Tzzz)

    snapShot2.Txxx[:] = fidMtx.MTxxxBTxxx .* snapShot1.Txxx + fidMtx.MTxxxBVx*(snapShot2.Vxx+snapShot2.Vxz)
    snapShot2.Txxz[:] = fidMtx.MTxxzBTxxz .* snapShot1.Txxz + fidMtx.MTxxzBVz*(snapShot2.Vzx+snapShot2.Vzz)
    snapShot2.Tzzx[:] = fidMtx.MTzzxBTzzx .* snapShot1.Tzzx + fidMtx.MTzzxBVx*(snapShot2.Vxx+snapShot2.Vxz)
    snapShot2.Tzzz[:] = fidMtx.MTzzzBTzzz .* snapShot1.Tzzz + fidMtx.MTzzzBVz*(snapShot2.Vzx+snapShot2.Vzz)
    snapShot2.Txzx[:] = fidMtx.MTxzxBTxzx .* snapShot1.Txzx + fidMtx.MTxzxBVz*(snapShot2.Vzx+snapShot2.Vzz)
    snapShot2.Txzz[:] = fidMtx.MTxzzBTxzz .* snapShot1.Txzz + fidMtx.MTxzzBVx*(snapShot2.Vxx+snapShot2.Vxz)
    return nothing
end

function OneStepAdjoint!(SnapShot2::SnapShot, SnapShot1::SnapShot, fidMtx::FidMtx)
    SnapShot2.it = SnapShot1.it - 1
    tmp = (fidMtx.MTxxxBVx)'*SnapShot1.Txxx + (fidMtx.MTzzxBVx)'*SnapShot1.Tzzx + (fidMtx.MTxzzBVx)'*SnapShot1.Txzz
    SnapShot2.Vxx = SnapShot1.Vxx + tmp
    SnapShot2.Vxz = SnapShot1.Vxz + tmp
    tmp = (fidMtx.MTxxzBVz)'*SnapShot1.Txxz + (fidMtx.MTzzzBVz)'*SnapShot1.Tzzz + (fidMtx.MTxzxBVz)'*SnapShot1.Txzx
    SnapShot2.Vzx = SnapShot1.Vzx + tmp
    SnapShot2.Vzz = SnapShot1.Vzz + tmp
    SnapShot2.Txxx = fidMtx.MTxxxBTxxx*SnapShot1.Txxx
    SnapShot2.Txxz = fidMtx.MTxxzBTxxz*SnapShot1.Txxz
    SnapShot2.Tzzx = fidMtx.MTzzxBTzzx*SnapShot1.Tzzx
    SnapShot2.Tzzz = fidMtx.MTzzzBTzzz*SnapShot1.Tzzz
    SnapShot2.Txzx = fidMtx.MTxzxBTxzx*SnapShot1.Txzx
    SnapShot2.Txzz = fidMtx.MTxzzBTxzz*SnapShot1.Txzz

    tmp = (fidMtx.MVxxBTxx)'*SnapShot2.Vxx
    SnapShot2.Txxx = tmp + SnapShot2.Txxx
    SnapShot2.Txxz = tmp + SnapShot2.Txxz
    tmp = (fidMtx.MVzzBTzz)'*SnapShot2.Vzz
    SnapShot2.Tzzx = tmp + SnapShot2.Tzzx
    SnapShot2.Tzzz = tmp + SnapShot2.Tzzz
    tmp = (fidMtx.MVxzBTxz)'*SnapShot2.Vxz + (fidMtx.MVzxBTxz)'*SnapShot2.Vzx
    SnapShot2.Txzx = tmp + SnapShot2.Txzx
    SnapShot2.Txzz = tmp + SnapShot2.Txzz
    SnapShot2.Vxx = fidMtx.MVxxBVxx*SnapShot2.Vxx
    SnapShot2.Vxz = fidMtx.MVxzBVxz*SnapShot2.Vxz
    SnapShot2.Vzx = fidMtx.MVzxBVzx*SnapShot2.Vzx
    SnapShot2.Vzz = fidMtx.MVzzBVzz*SnapShot2.Vzz
    return nothing
end

function OneStepAdjoint!(SnapShot2::SnapShot, SnapShot1::SnapShot, fidMtx::FidMtx)
    SnapShot2.it = SnapShot1.it - 1
    AcmulB!(tmp , fidMtx.MTxxxBVx, SnapShot1.Txxx); AcmulB!(tmp1, fidMtx.MTzzxBVx, SnapShot1.Tzzx); mysum2!(tmp, tmp1);
    AcmulB!(tmp1, fidMtx.MTxzzBVx, SnapShot1.Txzz); mysum2!(tmp , tmp1); mysum1!(SnapShot2.Vxx, SnapShot1.Vxx, tmp); mysum1!(SnapShot2.Vxz, SnapShot1.Vxz, tmp);

    AcmulB!(tmp , fidMtx.MTxxzBVz, SnapShot1.Txxz); AcmulB!(tmp1, fidMtx.MTzzzBVz, SnapShot1.Tzzz); mysum2!(tmp, tmp1);
    AcmulB!(tmp1, fidMtx.MTxzxBVz, SnapShot1.Txzx); mysum2!(tmp , tmp1); mysum1!(SnapShot2.Vzx, SnapShot1.Vzx, tmp); mysum1!(SnapShot2.Vzz, SnapShot1.Vzz, tmp);

    SnapShot2.Txxx = fidMtx.MTxxxBTxxx*SnapShot1.Txxx
    SnapShot2.Txxz = fidMtx.MTxxzBTxxz*SnapShot1.Txxz
    SnapShot2.Tzzx = fidMtx.MTzzxBTzzx*SnapShot1.Tzzx
    SnapShot2.Tzzz = fidMtx.MTzzzBTzzz*SnapShot1.Tzzz
    SnapShot2.Txzx = fidMtx.MTxzxBTxzx*SnapShot1.Txzx
    SnapShot2.Txzz = fidMtx.MTxzzBTxzz*SnapShot1.Txzz

    tmp = (fidMtx.MVxxBTxx)'*SnapShot2.Vxx
    SnapShot2.Txxx = tmp + SnapShot2.Txxx
    SnapShot2.Txxz = tmp + SnapShot2.Txxz
    tmp = (fidMtx.MVzzBTzz)'*SnapShot2.Vzz
    SnapShot2.Tzzx = tmp + SnapShot2.Tzzx
    SnapShot2.Tzzz = tmp + SnapShot2.Tzzz
    tmp = (fidMtx.MVxzBTxz)'*SnapShot2.Vxz + (fidMtx.MVzxBTxz)'*SnapShot2.Vzx
    SnapShot2.Txzx = tmp + SnapShot2.Txzx
    SnapShot2.Txzz = tmp + SnapShot2.Txzz
    SnapShot2.Vxx = fidMtx.MVxxBVxx*SnapShot2.Vxx
    SnapShot2.Vxz = fidMtx.MVxzBVxz*SnapShot2.Vxz
    SnapShot2.Vzx = fidMtx.MVzxBVzx*SnapShot2.Vzx
    SnapShot2.Vzz = fidMtx.MVzzBVzz*SnapShot2.Vzz
    return nothing
end
