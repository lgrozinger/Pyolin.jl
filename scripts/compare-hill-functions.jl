using Pyolin
using DataFrames
using CSV
using Plots
using StatsPlots
using LinearAlgebra
using Distributions
using LsqFit
using Symbolics
using LaTeXStrings

DATADIR = "/home/lewis/sauce/julia/Pyolin/data/"
FIGDIR = "/home/lewis/sauce/julia/Pyolin/figures/compare-hill-functions/"

px(mm) = mm * 96 / 25.4

frame = CSV.read(DATADIR * "experimentstats.csv", DataFrame)
groups = groupby(filter(r -> r.plasmid == "Lmra_n1", frame), [:strain, :backbone])

params = Matrix(undef, 4, length(groups))

plt = plot(
    guidefontsize=9,
    tickfontsize=9,
    legendfontsize=9,
    dpi=900,
    size=(px(100), px(60)),
    grid=false,
    xlabel=L"$R$",
    ylabel=L"$f_{YFP}(R)$",
)

for i in 1:length(groups)
    outputs   = Experiments(groups[i])
    inputs    = Inputs(outputs)
    standards = Standards(outputs)

    strain = outputs[1].strain
    backbone = outputs[1].backbone

    X = median.(rpuconvert(inputs, standards))
    Y = median.(rpuconvert(outputs, standards))
    h = Hill(rpuconvert(inputs, standards), rpuconvert(outputs, standards))
    params[:, i] .= [h.ymin, h.ymax, h.K, h.n]

    xs = range(0, 2.25, 256)
    ys = h.(xs)
    plot!(plt, xs, ys, label="$(strain) $(backbone)", linewidth=3)

    otherplot = notplot(rpuconvert(inputs, standards), rpuconvert(outputs, standards))
    savefig(otherplot, FIGDIR * "$(strain)-$(backbone)-$(outputs[1].plasmid)-notplot.svg")
end

savefig(plt, FIGDIR * "compare_hills.svg")
