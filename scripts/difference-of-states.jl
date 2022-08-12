using Pyolin
using DataFrames
using CSV
using Plots
using StatsPlots
using LinearAlgebra
using Distributions
using LaTeXStrings

DATADIR = "/home/lewis/sauce/julia/Pyolin/data/"
FIGDIR = "/home/lewis/sauce/julia/Pyolin/figures/statistical_functional/"

frame = CSV.read(DATADIR * "experimentstats.csv", DataFrame)
frame = filter(r -> r.plasmid ∉ [Pyolin.INPUT, Pyolin.STANDARD, Pyolin.AUTOFLUOR], frame)
frame = filter(r -> r.iptg == 0 || r.iptg == 2000, frame)
gates = groupby(frame, [:strain, :backbone, :plasmid])


for gate in gates
    e = Experiments(gate)
    diffs = sample(e[1], 8192) .- sample(e[2], 8192)
    failure = diffs[diffs .< 0]
    success = diffs[diffs .> 0]
    fn = "$(e[1].strain)-$(e[1].backbone)-$(e[1].plasmid)-diffs.svg"
    px(mm) = mm * 96 / 25.4
    plt = plot(size=(px(73), px(50)))
    histogram!(
        plt,
        failure,
        xlabel = L"YFP(0$\mu$M) - YFP(2000$\mu$M)",
        ylabel = "Frequency",
        guidefontsize=9,
        tickfontsize=9,
        legendfontsize=9,
        label = "$(round(length(failure) / length(diffs) * 100, digits=1))% Failure",
    )
    histogram!(
        plt,
        success,
        label = "$(round(length(success) / length(diffs) * 100, digits=1))% Success",
    )
    savefig(plt, FIGDIR * fn)

    gate.percent = repeat([length(success) / length(diffs)], 2)
end

results = unique(combine(gates, [:strain, :backbone, :plasmid, :percent]), [:strain, :backbone, :plasmid])

results.sixtysix = results.percent .> 0.66
results.eighty = results.percent .> 0.8
results.ninetyfive = results.percent .> 0.95
results.eightynine = results.percent .> 0.89

results.strain = convert.(String, results.strain)
results.backbone = convert.(String, results.backbone)
results.plasmid = convert.(String, results.plasmid)

contexts = groupby(results, [:strain, :backbone])
ninetyfive = [sum(context.ninetyfive) for context in contexts]
eightynine = [sum(context.eightynine) for context in contexts]
sixtysix = [sum(context.sixtysix) for context in contexts]

