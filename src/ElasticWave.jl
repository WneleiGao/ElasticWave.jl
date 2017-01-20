module ElasticWave
    using DSP, Requires
    const spmatveclib = abspath(joinpath(splitdir(Base.source_path())[1],"..","deps","builds","spmatvec"))
    include("dataStructure/dataStructure.jl")
    include("forwardAndAdjoint/forwardAndAdjoint.jl")
    include("spmatvec/spmatvec.jl")
    include("MOD/MOD.jl")
    @require PyPlot include("plotting/plotting.jl")
end
