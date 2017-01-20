using PyCall
@pyimport matplotlib.pyplot as plt
@pyimport matplotlib.lines as lines
@pyimport matplotlib.animation as anim

export SeisPlot,
       waveAnim

include("SeisPlot.jl")
include("waveAnim.jl")
