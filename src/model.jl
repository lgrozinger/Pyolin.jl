const UD = UnivariateDistribution
const MD = MultivariateDistribution

struct Standardisation
    strain::String
    backbone::String
    prior_k::UD
    prior_θ::UD
    e::Experiment
    chain::Chains
end
function Standardisation(e::Experiment, k::UD, θ::UD, N::Int=1024)
    y = Float64.(sample(e, N; replace=true))
    model = _standardisation(y, k, θ)
    Turing.setadbackend(:zygote)
    chain = sample(model, NUTS(N, 0.65), N * 8)
    Standardisation(e.strain, e.backbone, k, θ, e, chain)
end
Standardisations(es::Vector{Experiment}, k::UD, θ::UD, N::Int=1024) = Standardisation.(e, k, θ, N)

marginalk(s::Standardisation) = Array(group(s.chain, "k"), [:parameters])[:]
marginalθ(s::Standardisation) = Array(group(s.chain, "θ"), [:parameters])[:]
function StatsBase.sample(s::Standardisation, N::Int)
    chain = sample(s.chain, N)
    k = Array(group(chain, "k"), [:parameters])
    θ = Array(group(chain, "θ"), [:parameters])
    idxs = sample(1:length(k[:]), N; replace=false)
    rand.(Gamma.(k[:][idxs], θ[:][idxs]))
end

@model function _standardisation(y, kdist::UD, θdist::UD, N::Int=1)
    if y === missing
        y = Vector{Float64}(undef, N)
    end
    k ~ filldist(kdist, length(y))
    θ ~ filldist(θdist, length(y))
    y ~ arraydist(Gamma.(k, θ))
    return y
end

function plotit(s::Standardisation)
    opts = Dict(
        :guidefontsize => 10,
        :tickfontsize => 8,
        :legendfontsize => 10,
    )
    
    A = density(marginalk(s), xlabel="k", ylabel="Density", label="")
    B = density(marginalθ(s), xlabel=L"\theta", ylabel="Density", label="")
    C = histogram2d(
        marginalk(s),
        marginalθ(s),
        xlabel="k",
        ylabel=L"\theta",
        label="",
        cbar=false
    )
    D = histogram(sample(s.e, 4096); xlabel="y", ylabel="Density", label="Sampled", normalize=true)
    density!(D, sample(s, 4096), label="Estimated", linewidth=3)
    plot!(D, xlims=(0, quantile(s.e, 0.98)))

    l = @layout [
        [a b c]
        d
    ]
    plot(A, B, C, D, layout=l)
end

function hypothesistest(s::Standardisation)
    x = sample(s, 4096)
    y = sample(s.e, 4096)
    ApproximateTwoSampleKSTest(x, y)
end


@model function InputSensor(x::Vector{T}, y, μ::Vector{T}, σ::Matrix{T}) where {T<:Real}
    if y === missing
        y = Vector{T}(undef, length(x))
    end
    p ~ MvLogNormal(μ, σ)
    k = p[1] ./ (1 .+ (p[2] ./ x).^p[3]) .+ eps(T)
    y ~ arraydist(Gamma.(k, p[4]))
    return y
end
function InputSensor(es::Vector{<:Experiment}, N::Int, μ::Vector{T}, σ::Matrix{T}) where {T<:Real}
    x = Vector{T}(undef, length(es) * N)
    y = Vector{T}(undef, length(es) * N)
    for i in 1:length(es)
        e = es[i]
        x[N * i - N + 1 : N * i] .= e.iptg
        y[N * i - N + 1 : N * i] .= sample(events(e)[CHANNEL], N; replace=false)
    end
    InputSensor(x, y, μ, σ)
end
function InputSensor(c::Chains)
    p = Array(c, [:parameters])
    p = fit(MvNormal, log.(p'))
    x -> InputSensor(x, missing, p.μ, collect(p.Σ))
end

@model function Inverter(i::Vector{T}, y, μ::Vector{T}, σ::Matrix{T}) where {T<:Real}
    if y === missing
        y = Vector{T}(undef, length(i))
    end
    p ~ MvLogNormal(μ, σ)
    a = p[1] ./ (1 .+ (p[2] ./ i).^p[3]) .+ eps(T)
    x ~ arraydist(Gamma.(a, p[4]))
    b = p[5] ./ (1 .+ (x ./ p[6]).^p[7]) .+ eps(T)
    y ~ arraydist(Gamma.(b, p[8]))
    return x, y
end
function Inverter(es::Vector{<:Experiment}, N::Int, μ::Vector{T}, σ::Vector{T}) where {T<:Real}
    x = Vector{T}(undef, length(es) * N)
    y = Vector{T}(undef, length(es) * N)
    for i in 1:length(es)
        e = es[i]
        x[N * i - N + 1 : N * i] .= e.iptg
        y[N * i - N + 1 : N * i] .= sample(events(e)[CHANNEL], N; replace=false)
    end
    Inverter(x, y, μ, σ)
end
function Inverter(c::Chains)
    p = Array(c, [:parameters])
    p = fit(MvNormal, log.(p'))
    i -> Inverter(i, missing, p.μ, p.Σ)
end
