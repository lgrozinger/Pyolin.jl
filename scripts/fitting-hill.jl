using Pyolin
using DataFrames
using CSV
using Plots
using StatsPlots
using LinearAlgebra
using Distributions
using LsqFit
using LaTeXStrings
using Symbolics

DATADIR = "/home/lewis/sauce/julia/Pyolin/data/"
FIGDIR = "/home/lewis/sauce/julia/Pyolin/figures/fitting-hill/"

frame = CSV.read(DATADIR * "experimentstats.csv", DataFrame)


xs = filter(r -> r.plasmid == "1818", frame)
xs = filter(r -> r.strain == "KT2440", xs)
xs = filter(r -> r.backbone == "pSeva221", xs)
sort!(xs, [:iptg])

ys = filter(r -> r.plasmid == "Lmra_n1", frame)
ys = filter(r -> r.strain == "KT2440", ys)
ys = filter(r -> r.backbone == "pSeva221", ys)
sort!(ys, [:iptg])

inputs = Experiments(xs)
outputs = Experiments(ys)
autos = Autos(inputs)
standards = Standards(inputs)

px(mm) = mm * 96 / 25.4

opts = Dict(
    :guidefontsize => 9,
    :tickfontsize => 9,
    :legendfontsize => 4,
    :legend => :topright,
    :grid => false,
    :dpi => 900,
)

rpuinputs = rpuevents(inputs, standards)
rpuoutputs = rpuevents(outputs, standards)

pltA = plot(
    ;xlabel="YFP",
    ylabel="Density",
    xlims=(0, 1.5),
    size=(px(70), px(42)),
    opts...,
)

for (i, e) in enumerate(inputs)
    density!(pltA, rpuinputs[i], label="$(e.iptg)" * L"$\mu$M", linewidth=2, color=i)
    vline!(pltA, [median(rpuinputs[i])], color=i, label="", linestyle=:dash)
end
savefig(pltA, FIGDIR * "fitting-hill-A.svg")

pltB = plot(
    ;xlabel="YFP",
    ylabel="Density",
    xlims=(0, 4),
    size=(px(73), px(43)),
    opts...,
)

for (i, e) in enumerate(outputs)
    density!(pltB, rpuoutputs[i], label="$(e.iptg)" * L"$\mu$M", linewidth=2, color=i)
    vline!(pltB, [median(rpuoutputs[i])], color=i, label="", linestyle=:dash)
end
savefig(pltB, FIGDIR * "fitting-hill-B.svg")

iptgs = getproperty.(inputs, :iptg)
X = median.(rpuconvert(inputs, standards))
Y = median.(rpuconvert(outputs, standards))

pltC = plot()
plot!(pltC, iptgs, X, color=1, label=false, linewidth=3)
plot!(pltC, iptgs, Y, color=2, label=false, linewidth=3)
scatter!(
    pltC,
    iptgs,
    X,
    xlabel="IPTG",
    ylabel="YFP",
    label=L"$R_i$",
    legendfontsize=10,
    guidefontsize=10,
    markersize=5,
    dpi=900,
    linewidth=2,
    size=(px(73), px(43)),
    grid=false,
    color=1
)
scatter!(
    pltC,
    iptgs,
    Y,
    xlabel="IPTG",
    ylabel="YFP",
    label=L"$f_{YFP}$",
    legendfontsize=10,
    guidefontsize=10,
    markersize=5,
    dpi=900,
    linewidth=2,
    size=(px(73), px(43)),
    grid=false,
    color=2
)

savefig(pltC, FIGDIR * "fitting-hill-C.svg")

h = Hill(rpuconvert(inputs, standards), rpuconvert(outputs, standards))
pltD = plot(
    xlabel=L"$R$",
    ylabel=L"$f_{YFP}(R)$",
    guidefontsize=10,
    legendfontsize=10,
    dpi=900,
    linewidth=2,
    size=(px(73), px(43)),
    grid=false,
)

xs = range(0, 1, 100)
ys = h.(xs)
plot!(pltD, xs, ys, label="Fitted model", linewidth=3)
scatter!(pltD, X, Y, label="Data points", markersize=3)
savefig(pltD, FIGDIR * "fitting-hill-D.svg")
