using Pyolin
using DataFrames
using CSV
using LaTeXStrings
using Plots
using StatsPlots
using LinearAlgebra
using Distributions
using InlineStrings


DATADIR = "data/"
FIGDIR = "figures/chains/"

frame = CSV.read(pwd()*"/"*DATADIR*"experimentstats.csv", DataFrame)
frame = filter!(r -> r.plasmid in ["1818"], frame)
E = Experiments(frame)

function work(e::Experiment)
    C = Constitutive(e, 256)
    save(C, pwd()*"/"*DATADIR*"chains/$(e.strain)-$(e.backbone)-$(e.plasmid)-$(e.iptg)-const.h5")
    plt = constitutiveplot(C)
    savefig(plt, pwd()*"/"*FIGDIR*"$(e.strain)-$(e.backbone)-$(e.plasmid)-$(e.iptg).svg")
    g = fit(Gamma, e)
    p = posteriorpredict(C, 2048)
    [
        e.strain, e.backbone, e.plasmid, e.iptg,
        mean(e), var(e),
        mean(simulate(C, 2048)), var(simulate(C, 2048)),
        mean(g), var(g),
        mean(p[1, :]), mean(p[2, :]),
        var(p[1, :]), var(p[2, :]),
    ]
end

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
        meanα=Float32[],
        meanβ=Float32[],
        varα=Float32[],
        varβ=Float32[],
    )

    for (i, e) in enumerate(E)
        println("Iteration $i")
        append!(results, DataFrame([[r] for r in work(e)], propertynames(results)))
    end
    results
end

results = process()
results |> CSV.write(pwd()*"/"*DATADIR*"chains/constitutive-chains.csv")


