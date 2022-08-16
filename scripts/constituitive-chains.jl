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


DATADIR = "data/"
FIGDIR = "figures/chains/"

frame = CSV.read(pwd()*"/"*DATADIR*"experimentstats.csv", DataFrame)
frame = filter!(r -> r.plasmid != "1201", frame)
E = Experiments(frame)

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

    function work(e)
        tries = 0
        C = nothing
        while tries < 16
            try
                C = Constitutive(e, 512)
                tries = 9999
            catch DomainError
                tries = tries + 1
            end
        end
        save(C, pwd()*"/"*DATADIR*"chains/$(e.strain)-$(e.backbone)-$(e.plasmid)-$(e.iptg)-const.h5")
        plt = constitutiveplot(C)
        savefig(plt, pwd()*"/"*FIGDIR*"$(e.strain)-$(e.backbone)-$(e.plasmid)-$(e.iptg).svg")
        g = fit(Gamma, e)
        p = posteriorpredict(C, 2048)
        row = [
            e.strain, e.backbone, e.plasmid, e.iptg,
            mean(e), var(e),
            mean(simulate(C, 2048)), var(simulate(C, 2048)),
            mean(g), var(g),
            mean(p[1, :]), mean(p[2, :]),
            var(p[1, :]), var(p[2, :]),
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
results |> CSV.write(pwd()*"/"*DATADIR*"chains/constitutive-chains.csv")


