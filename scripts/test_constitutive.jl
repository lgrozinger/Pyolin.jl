using Pyolin
using DataFrames
using CSV
using Plots
using StatsPlots
using LinearAlgebra
using Distributions


DATADIR = "/home/lewis/sauce/julia/Pyolin/data/"
FIGDIR = "/home/lewis/sauce/julia/Pyolin/figures/dispersion/"

frame = CSV.read(DATADIR * "experimentstats.csv", DataFrame)
frame = filter(r -> r.plasmid == "1717", frame)

E = Experiments(frame)





