using Pyolin
using DataFrames
using ProgressMeter
using CSV
using LaTeXStrings
using Plots
using StatsPlots
using LinearAlgebra
using Distributions
using InlineStrings


DATADIR = "/home/lewis/sauce/julia/Pyolin/data/"
FIGDIR = "/home/lewis/sauce/julia/Pyolin/figures/dispersion/"

frame = CSV.read(DATADIR * "experimentstats.csv", DataFrame)
E = Experiments(frame[103:104, :])

function process()
    results = DataFrame(
        strain=String63[],
        backbone=String63[],
        plasmid=String63[],
        iptg=Int[],
        samplemean=Float32[],
        samplevar=Float32[],
        postmean=Float32[],
        postvar=Float32[],
        exactmean=Float32[],
        exactvar=Float32[],
        estimvar=Float32[]
    )

    function work(e)
        C = Constitutive(e, 512)
        save(C, DATADIR * "chains/$(e.strain)-$(e.backbone)-$(e.plasmid)-$(e.iptg)-const.h5")
        plt = constitutiveplot(s)
        savefig(plt, FIGDIR * "$(e.strain)-$(e.backbone)-$(e.plasmid)-$(e.iptg).svg")
        g = fit(Gamma, e)
        row = [
            e.strain, e.backbone, e.plasmid, e.iptg,
            mean(e), var(e),
            mean(simulate(C, 2048)), var(simulate(C, 2048)),
            mean(g), var(g)
        ]
        DataFrame([[r] for r in row], propertynames(results))
    end

    for (i, e) in enumerate(E)
        println("Iteration $i")
        append!(results, work(e))
    end
    results
end

results = process()

