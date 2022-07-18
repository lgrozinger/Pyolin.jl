using Pyolin
using Turing
using MCMCChains
using Plots
using StatsPlots
using Plots.Measures
using DataFrames
using CSV
using LinearAlgebra
using ProgressMeter
using Random
using LaTeXStrings
using HypothesisTests


const DATADIR = "/home/campus.ncl.ac.uk/b8051106/sauce/julia/Pyolin.jl/data/"
const FIGDIR = "/home/campus.ncl.ac.uk/b8051106/sauce/julia/Pyolin.jl/figures/"

commonopts = Dict(
    :guidefontsize => 12,
    :legend => false,
    :tickfontsize => 12,
    :markersize => 6,
    :markeralpha => 0.5,
    :markerstrokealpha => 0.5
)

function densitysubplot(data, letter; kwargs...)
    opts = Dict(
        :ylabel => "Density",
        :label => false,
        :linewidth => 2,
        :tickfontsize => 8,
        :titlefontsize => 10,
        :titleloc => :left,
        :title => letter
    )
    density(data; merge(opts, kwargs)...)
end

function plotstandardisationchain(chain)
    f = Standardisation(chain)
    p = Array(chain, [:parameters])
    Pp = fit(MvNormal, log.(p'))
    Pp = MvLogNormal(Pp.μ, collect(Pp.Σ))

    krange = range(minimum(p[:,1]), maximum(p[:,1]), 1024)
    θrange = range(minimum(p[:,2]), maximum(p[:,2]), 1024)
    kθ = heatmap(
        krange, θrange, Pp;
        titleloc=:left,
        title="C",
        titlefontsize=10,
        cbar=false,
        tickfontsize=8,
        ylabel=L"\theta",
        xlabel="k"
    )
    k = densitysubplot(p[:,1], "A"; xlabel="k")
    θ = densitysubplot(p[:,2], "B"; xlabel=L"\theta", ylabel="")
    s = reduce(vcat, f() for _ in 1:2^16)
    y = densitysubplot(s, "D"; xlabel="y", label="Estimated", tickfontsize=10, linewidth=3)
    l = @layout [
        [a b c]
        d
    ]
    plot(k, θ, kθ, y, layout=l)
end

function plotinputchain2d(chain; kwargs...)
    f = InputSensor(chain)
    x = eps(Float64) .+ 2000 .* rand(2^14)
    y = f(x)()
    histogram2d(x, y; kwargs...)
end

function plotinputchain(chain, x; kwargs...)
    f = InputSensor(chain)
    y = f(repeat([x], 2^14))()
    density(y; kwargs...)
end

frame = CSV.read(DATADIR * "experimentstats.csv", DataFrame)
standardisations = filter(r -> r.plasmid == "1717", frame)
inputsensors = filter(r -> r.plasmid == "1818" && r.iptg != 0, frame)
