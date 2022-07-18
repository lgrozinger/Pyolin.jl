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

function standardisationchain(experiment; N=8192, M=8192)
    μ = [1.0, 1.0]
    Σ = Matrix(0.1I, 2, 2)
    model = Standardisation(experiment, N, μ, Σ)
    Turing.setadbackend(:forwarddiff)
    sample(model, NUTS(M, 0.8), M)
end
function inputchain(experiments)
    μ = [1.0, 1.0, 1.0, 1.0]
    Σ = Matrix(0.1I, 4, 4)
    model = InputSensor(experiments, 512, μ, Σ)
    Turing.setadbackend(:forwarddiff)
    sample(model, NUTS(512, 0.65), 4096)
end
function inverterchain(experiments)
    μ = repeat([1.0], 8)
    Σ = Matrix(0.1I, 8, 8)
    model = Inverter(experiments, 256, μ, Σ)
    Turing.setadbackend(:zygote)
    sample(model, NUTS(512, 0.65), 4096)
end

function contouring(x, y)
    f = fit(MvNormal, log.(hcat(x, y)))
    f = MvLogNormal(f.μ, collect(f.Σ))
    ex = minimum(x):0.01:maximum(x)
    why = minimum(y):0.01:maximum(y)
    X = repeat(reshape(ex, 1, :), length(y), 1)
    Y = repeat(y, 1, length(x))
    Z = map((x, y) -> pdf(f, [x, y]), X, Y)
    contour(x, y, Z; cbar=false)
end

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

function __savechain__(experiment, chain)
    fn = "chain-$(experiment.strain)-$(experiment.backbone)-$(experiment.plasmid)-$(experiment.iptg).csv"
    DataFrame(chain) |> CSV.write(fn)
end

frame = CSV.read(DATADIR * "experimentstats.csv", DataFrame)
standardisations = filter(r -> r.plasmid == "1717", frame)
inputsensors = filter(r -> r.plasmid == "1818" && r.iptg != 0, frame)

function process_standardisation(e)
    chain = standardisationchain(e)
    __savechain__(e, chain)
    plt = plotstandardisationchain(chain)
    samples = sample(events(e)[Pyolin.CHANNEL], 2^15)
    histogram!(plt, samples; label="Sampled", subplot=4, normalize=:true)

    f = Standardisation(chain)
    fsamples = reduce(vcat, f() for _ in 1:2^15)

    ks = ApproximateTwoSampleKSTest(samples, fsamples)
    turingstat = sqrt(ks.n_x * ks.n_y / (ks.n_x + ks.n_y)) * ks.δ
    turingpvalue = pvalue(ks)
    @show turingstat
    @show turingpvalue
    ks = ApproximateTwoSampleKSTest(samples, rand(fit(Gamma, e), 2^15))
    fitstat = sqrt(ks.n_x * ks.n_y / (ks.n_x + ks.n_y)) * ks.δ
    fitpvalue = pvalue(ks)
    @show fitstat
    @show fitpvalue

    savefig(plt, FIGDIR * "standardisations/$(e.strain)-$(e.backbone)-$(e.iptg).svg")
    turingstat, turingpvalue, fitstat, fitpvalue
end

function process_standardisations(frame)
    scores = map(process_standardisation, Experiments(frame))
    scores
end
