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
standardisations = filter(r -> r.plasmid == "1717", frame)
inputsensors = filter(r -> r.plasmid == "1818", frame)
nots = filter(r -> r.plasmid ∉ ["1818", "1717", "1201"], frame)

E = Experiments(vcat(standardisations, inputsensors, nots))
prior = MvLogNormal(ones(2), Matrix(0.5I, 2, 2))

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
        s = Constituitive(e, prior, 2048)
        savechain(s, DATADIR * "chains/")
        
        plt = constituitiveplot(s)
        savefig(plt, FIGDIR * "$(strain(s))-$(backbone(s))-$(plasmid(s))-$(iptg(s))" * ".svg")
        g = fit(Gamma, e)
        row = [
            strain(s), backbone(s), plasmid(s), iptg(s),
            mean(e), var(e),
            mean(s), var(s),
            mean(g), var(g),
            varvar(s, 2048)
        ]
        DataFrame([[r] for r in row], propertynames(results))
    end

    for (i, e) in enumerate(E)
        println("Iteration $i")
        append!(results, work(e))
    end
    results
end

# process()

contexts = groupby(standardisations, [:strain, :backbone])
px(mm) = mm * 96 / 25.4
splt = plot(size=(px(100), px(60)), yscale=:log10)
for (i, context) in enumerate(contexts[1:3])
    x = [Constituitive(e, prior, DATADIR * "chains/") for e in Experiments(context)]
    constituitivecov!(splt, x; n_std=2.5, seriescolor=i)
end

# s = Constituitive(Experiments(contexts[1])[1], prior, 2048)

# hist = plot(size=(px(100), px(60)))
# histogram2d!(hist, ksamples(s), θsamples(s); colorbar=false, xlabel="k", ylabel=L"\theta")

# example = plot(size=(px(200), px(60)))
# histogram!(example, sample(experiment(s), 4096), label="Sampled", ylabel="Density", normalize=true)
# density!(example, sample(s, 4096), label="Estimated", xlabel="YFP", ylabel="Density", linewidth=3)

contexts = groupby(inputsensors, [:strain, :backbone])
px(mm) = mm * 96 / 25.4
iplt = plot(size=(px(100), px(60)))
for (i, context) in enumerate(contexts[1:3])
    x = [Constituitive(e, prior, DATADIR * "chains/") for e in Experiments(context)]
    constituitivecov!(iplt, x; n_std=2.5, seriescolor=i, yscale=:log10, xscale=:log10)
end

contexts = groupby(nots, [:strain, :backbone, :plasmid])
px(mm) = mm * 96 / 25.4
nplt = plot(size=(px(100), px(60)))
for (i, context) in enumerate(contexts[1:3])
    x = [Constituitive(e, prior, DATADIR * "chains/") for e in Experiments(context)]
    constituitivecov!(nplt, x; n_std=2.5, seriescolor=i, yscale=:log10)
end
