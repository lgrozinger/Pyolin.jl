using Pyolin
using Turing
using MCMCChains
using Plots
using Plots
using Plots.Measures
using DataFrames
using CSV
using LinearAlgebra
using ProgressMeter
using Random


const DATADIR = "/home/lewis/sauce/julia/Pyolin/data/"
const FIGDIR = "/home/lewis/sauce/julia/Pyolin/figures/"

commonopts = Dict(
    :guidefontsize => 12,
    :legend => false,
    :tickfontsize => 12,
    :markersize => 6,
    :markeralpha => 0.5,
    :markerstrokealpha => 0.5
)

function __statisticpoint__(xmodel, ymodel, statsfun, N)
    xsamples = reduce(vcat, xmodel() for _ in 1:N)
    ysamples = reduce(vcat, ymodel() for _ in 1:N)
    statsfun(xsamples), statsfun(ysamples)
end

function __statisticpoints__(xmodels, ymodels, statsfun, N)
    [__statisticpoint__(x, y, statsfun, N) for (x, y) in zip(xmodels, ymodels)]
end

function __statisticplot__(xmodels, ymodels, statsfun, N; kwargs...)
    xy = __statisticpoints__(xmodels, ymodels, statsfun, N)
    plt = scatter(first.(xy), last.(xy); kwargs...)
    min, max = minimum(first.(xy)), maximum(first.(xy))
    plot!(plt, [min, max], [min, max], label=false, linestyle=:dash)
    plt
end

function __checktypes__(model)
    @code_warntype model.f(
        model,
        Turing.VarInfo(model),
        Turing.SamplingContext(
            Random.GLOBAL_RNG, Turing.SampleFromPrior(), Turing.DefaultContext(),
        ),
        model.args...,
    )
end

function __sampleexperiment__(experiment)
    evs = events(experiment)
    () -> sample(evs[Pyolin.CHANNEL], 1)
end

function __sampledistribution__(experiment; type=Gamma)
    distribution = fit(type, experiment)
    () -> rand(distribution, 1)
end

function __savechain__(experiment, chain)
    fn = "chain-$(experiment.strain)-$(experiment.backbone)-$(experiment.plasmid)-$(experiment.iptg).csv"
    DataFrame(chain) |> CSV.write(fn)
end

function __samplestandard__(experiment, N=4096; witherror=true)
    priors = ([1.0, 1.0], Matrix(0.25I, 2, 2))
    model = model1717(experiment, N, priors...; witherror=witherror)
    chain = sample(model, NUTS(256, 0.65), 4096; drop_warmup=true, verbose=false)
    __savechain__(experiment, chain)
    psamples = reduce(hcat, x[:] for x in get_params(chain).p)
    mvn = fit(MvNormal, log.(psamples'))
    if witherror
        model1717(1, mvn.μ, collect(mvn.Σ), model.args.ϵ)
    else
        model1717(1, mvn.μ, collect(mvn.Σ))
    end
end

frame = CSV.read(DATADIR * "experimentstats.csv", DataFrame)
standardisations = filter(r -> r.plasmid == "1717", frame)
inputsensors = filter(r -> r.plasmid == "1818" && r.iptg != 0, frame)

function __bayesstatistics__(experiments)
    bayes = [__samplestandard__(experiment; witherror=false) for experiment in experiments]
    experimental = [__sampleexperiment__(experiment) for experiment in experiments]

    means = __statisticpoints__(bayes, experimental, mean, 2^14)
    variances = __statisticpoints__(bayes, experimental, var, 2^14)
    means, variances
end
# mom_mean(Experiments(standardisations), "means-1717-method-of-moments.svg")
# mom_var(Experiments(standardisations), "vars-1717-method-of-moments.svg")

# mom_mean(Experiments(inputsensors), "means-1818-method-of-moments.svg")
# mom_var(Experiments(inputsensors), "vars-1818-method-of-moments.svg")

function bayesian_standard(experiments, outputfn)
    N = 512
    M = 4096

    Ne = length(experiments)
    x = Matrix{Float64}(undef, 2, Ne)
    y = Matrix{Float64}(undef, 2, Ne)
    prior = ([1.0, 1.0], Matrix(0.1I, 2, 2))

    for i in 1:Ne
        samples = sample(events(experiments[i])[Pyolin.CHANNEL], N; replace=false)
        model = StandardModel(model1717(samples, prior...), M)
        y[1, i] = mean(samples)
        y[2, i] = var(samples)
        posterior = model(N)
        x[1, i] = mean(posterior[3, :])
        x[2, i] = var(posterior[3, :])
    end
        
    opts = Dict(
        :xlabel => "Model mean",
        :ylabel => "Sample mean",
        :guidefontsize => 12,
        :legend => false,
        :tickfontsize => 12,
        :markersize => 6
    )
    plt = scatter(x[1, :], y[1, :]; opts...)
    line = [minimum(x[1, :]), maximum(x[1, :])]
    plot!(plt, line, line)
    savefig(plt, FIGDIR * "models/means-" * outputfn)

    opts = merge(
        opts,
        Dict(:xlabel=>"Model variance", :ylabel=>"Sample variance")
    )
    plt = scatter(x[2, :], y[2, :]; opts...)
    line = [minimum(x[2, :]), maximum(x[2, :])]
    plot!(plt, line, line, xscale=:log10, yscale=:log10)
    savefig(plt, FIGDIR * "models/vars-" * outputfn)
end

function bayesian_inputsensor(experiments, outputfn)
    N = 512
    M = 4096
    prior = ([1.0, 1.0, 1.0, 1.0], Matrix(0.5I, 4, 4))
    model = cond1818(experiments, N, prior...)
    InputModel(model, M)
end

# bayesian_standard(Experiments(standardisations), "1717-bayesian-model.svg")
# m = bayesian_inputsensor(Experiments(filter(r -> r.strain == "KT2440" && r.backbone == "pSeva221", inputsensors)), nothing)
