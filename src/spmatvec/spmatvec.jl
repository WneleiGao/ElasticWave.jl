export ParAmulB!,
       AmulB!,
       AcmulB!,
       ca_add_b!

include("AmulB!.jl")
include("AcmulB!.jl")
include("a_add_b.jl")


m = 1000
A = sprand(m, m, 0.9);
x = rand(m);
y = zeros(m);
y1 = zeros(m);
nthreads = 4;
yt = zeros(m*nthreads);

@time ParAmulB!(nthreads, yt, y, A, x);
@time A_mul_B!(y1, A, x);
