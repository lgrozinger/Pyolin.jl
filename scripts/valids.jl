using Pyolin
using DataFrames
using CSV


DATADIR = "/home/lewis/sauce/julia/Pyolin/data/"
FIGDIR  = "/home/lewis/sauce/julia/Pyolin/figures/"

frame = CSV.read(DATADIR * "gatestatsrpu.csv", DataFrame)
contexts = groupby(frame, [:strain, :backbone])

for context in contexts 
    gates = Hills(context)
    println("$(context.strain[1]), $(context.backbone[1]): ", sum(Pyolin.valid.(gates)))
end


