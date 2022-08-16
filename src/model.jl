abstract type Model end

# to implement
function simulate end
function simulateprior end
function priorpredict end
function posteriorpredict end
function conditions end

chain(x::Model) = x.chain
model(x::Model) = x.model

function save(x::Model, fn::AbstractString)
    h5open(fn, "w") do f
        write(f, chain(x))
    end
end

# a constitutive expression model
@model function constitutive(μ::LogNormal{T}, σ::LogNormal{T}, y::Vector{T}) where {T<:Real}
    μ₁ ~ μ
    μ₂ ~ μ
    σ₁ ~ σ
    σ₂ ~ σ
    α ~ filldist(LogNormal(μ₁, σ₁), length(y))
    β ~ filldist(LogNormal(μ₂, σ₂), length(y))
    y ~ arraydist(Gamma.(α, β))
    return y
end
@model function constitutive(μ::LogNormal{T}, σ::LogNormal{T}, N::Int) where {T<:Real}
    y = Vector{T}(undef, N)
    μ₁ ~ μ
    μ₂ ~ μ
    σ₁ ~ σ
    σ₂ ~ σ
    α ~ filldist(LogNormal(μ₁, σ₁), length(y))
    β ~ filldist(LogNormal(μ₂, σ₂), length(y))
    y ~ arraydist(Gamma.(α, β))
end

struct Constitutive <: Model
    chain::Chains
    model::DynamicPPL.Model
    e    ::Experiment
end

function Constitutive(e::Experiment, μ::LogNormal{T}, σ::LogNormal{T}, N::Int=2048) where {T<:Real}
    y = T.(sample(e, N; replace=false))
    model = constitutive(μ, σ, y)
    Turing.setadbackend(:zygote)
    chain = sample(model, DynamicNUTS(), N * 8, save_state=true, discard_adapt=true)
    Constitutive(chain, model, e)
end

function Constitutive(e::Experiment, N::Int=2048)
    Constitutive(e, LogNormal(2f0, 5f-1), LogNormal(2f0, 5f-1), N)
end

function Constitutive(e::Experiment, fn::AbstractString)
    chain = h5open(fn, "r") do f
        read(f, Chains)
    end
    Constitutive(chain, chain.info.model, e)
end

function priorpredict(x::Constitutive, N::Int)
    hcat(
        rand.(LogNormal.(rand(model(x).args.μ, N), rand(model(x).args.σ, N))),
        rand.(LogNormal.(rand(model(x).args.μ, N), rand(model(x).args.σ, N))),
    )'
end

function priorsimulate(x::Constitutive, N::Int)
    p = priorpredict(x, N)
    rand.(Gamma.(p[1, :], p[2, :]))
end

function posteriorpredict(x::Constitutive, N::Int)
    n = length(chain(x))
    idxs = sample(1:n, N; replace=true)
    μ₁ = Array(group(chain(x), "μ₁"), [:parameters])[:][idxs]
    μ₂ = Array(group(chain(x), "μ₂"), [:parameters])[:][idxs]
    σ₁ = Array(group(chain(x), "σ₁"), [:parameters])[:][idxs]
    σ₂ = Array(group(chain(x), "σ₂"), [:parameters])[:][idxs]
    hcat(rand.(LogNormal.(μ₁, σ₁)), rand.(LogNormal.(μ₂, σ₂)))'
end

function simulate(x::Constitutive, N::Int)
    p = posteriorpredict(x, N)
    rand.(Gamma.(p[1, :], p[2, :]))
end

function load(fn::AbstractString, experiment::Experiment)
    chain = h5open(fn, "r") do f
        read(f, Chains)
    end
    model = chain.info.model
    Constitutive(chain, model, experiment)
end

struct Induced <: Model
    chain::Chains
    model::DynamicPPL.Model
    e    ::Vector{<:Experiment}
end

@userplot ConstitutivePlot
@recipe function f(x::ConstitutivePlot)
    c, = x.args

    guidefontsize --> 9
    tickfontsize --> 9
    legendfontsize --> 9
    titlefontsize --> 10

    l = @layout [
        [a b c]
        d
    ]
    layout --> l

    ps = posteriorpredict(c, 4096)

    @series begin
        seriestype := :density
        xlabel := L"$\alpha$"
        ylabel := "Density"
        label := false
        subplot := 2
        title := ""
        ps[1, :]
    end

    @series begin
        seriestype := :density
        xlabel := L"$\beta$"
        ylabel := "Density"
        label := false
        subplot := 3
        title := ""
        ps[2, :]
    end

    @series begin
        seriestype := :histogram2d
        bins := 100
        xlabel := L"$\alpha$"
        ylabel := L"$\beta$"
        subplot := 1
        colorbar := false
        title := ""
        ps[1,:], ps[2,:]
    end

    @series begin
        seriestype := :histogram
        normalize := true
        xlabel := "y"
        ylabel := "Density"
        label := "Sampled"
        subplot := 4
        sample(c.e, 4096)
    end

    @series begin
        seriestype := :density
        label := "Estimated"
        xlabel := "YFP"
        ylabel := "Density"
        linewidth := 3
        subplot := 4
        simulate(c, 8192)
    end
end


