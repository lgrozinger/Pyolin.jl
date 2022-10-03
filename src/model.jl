abstract type Model end

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
@model function constitutive(y::Vector{T}) where {T<:Real}
    # parameters
    μ₁ ~ LogNormal()
    μ₂ ~ LogNormal()
    μ₃ ~ LogNormal()
    σ₁ ~ LogNormal()
    σ₂ ~ LogNormal()
    σ₃ ~ LogNormal()

    # variables
    α ~ filldist(LogNormal(μ₁, σ₁), length(y))
    β ~ filldist(LogNormal(μ₂, σ₂), length(y))
    e ~ filldist(LogNormal(μ₃, σ₃), length(y))
    x ~ arraydist(Gamma.(α, β))
    y ~ arraydist(Normal.(x .+ e, 0.01))
    return y
end

struct Constitutive <: Model
    chain::Chains
    model::DynamicPPL.Model
    e    ::Experiment
end

function Constitutive(e::Experiment, N::Int=2048)
    y = Float64.(sample(e, N; replace=false))
    model = constitutive(y)
    Turing.setadbackend(:reversediff)
    Turing.setrdcache(true)
    opts = Dict(
        :save_state => true,
        :discard_adapt => true,
    )
    chain = sample(model, HMC(1e-6, 10), N * 32; opts...)
    Constitutive(chain, model, e)
end

function Constitutive(e::Experiment, fn::AbstractString)
    chain = h5open(fn, "r") do f
        read(f, Chains)
    end
    Constitutive(chain, chain.info.model, e)
end

function priorpredict(x::Constitutive, N::Int)
    hcat(
        rand.(LogNormal.(rand(LogNormal(), N)), rand(LogNormal(), N)),
        rand.(LogNormal.(rand(LogNormal(), N)), rand(LogNormal(), N)),
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
    μ₃ = Array(group(chain(x), "μ₃"), [:parameters])[:][idxs]
    σ₁ = Array(group(chain(x), "σ₁"), [:parameters])[:][idxs]
    σ₂ = Array(group(chain(x), "σ₂"), [:parameters])[:][idxs]
    σ₃ = Array(group(chain(x), "σ₃"), [:parameters])[:][idxs]
    hcat(rand.(LogNormal.(μ₁, σ₁)), rand.(LogNormal.(μ₂, σ₂)), rand.(LogNormal.(μ₃, σ₃)))'
end

function simulate(x::Constitutive, N::Int)
    p = posteriorpredict(x, N)
    rand.(Gamma.(p[1, :], p[2, :])) .+ p[3, :]
end


# an induced expression model
@model function induced(p::LogNormal{T}, x::Vector{T}, y::Vector{T}) where {T<:Real}
    β₁ ~ filldist(p, length(y))
    k₁ ~ filldist(p, length(y))
    ϵ₁ ~ filldist(p, length(y))
    n₁ ~ filldist(p, length(y))
    β ~ filldist(p, length(y))
    
    α = ϵ₁ .+ (β₁ ./ (1 .+ (x ./ k₁).^n₁))
    y ~ arraydist(Gamma.(α, β))
    return y
end
@model function induced(p::LogNormal{T}, x::Vector{T}, N::Int) where {T<:Real}
    y = Vector{T}(undef, N)
    β₁ ~ p
    k₁ ~ p
    ϵ₁ ~ p
    n₁ ~ p
    β ~ filldist(p, length(y))
    
    α = ϵ₁ .+ (β₁ ./ (1 .+ (x ./ k₁).^n₁))
    y ~ arraydist(Gamma.(α, β))
    return y
end

struct Induced <: Model
    chain::Chains
    model::DynamicPPL.Model
    e    ::Vector{<:Experiment}
end

function Induced(es::Vector{<:Experiment}, p::LogNormal{T}, N::Int=256) where {T<:Real}
    x = repeat(getproperty.(es, :iptg), inner=N)
    y = reduce(vcat, T.(sample(e, N; replace=false)) for e in es)
    model = induced(p, x, y)
    Turing.setadbackend(:zygote)
    chain = sample(model, DynamicNUTS(), N * 8, save_state=true, discard_adapt=true)
    Induced(chain, model, es)
end

function Induced(es::Vector{<:Experiment}, N::Int=256)
    Induced(es, LogNormal(2f0, 5f-1), N)
end

function Induced(es::Vector{<:Experiment}, fn::AbstractString)
    chain = h5open(fn, "r") do f
        read(f, Chains)
    end
    Induced(chain, chain.info.model, es)
end

function priorpredict(x::Induced, N::Int)
    hcat(
        rand.(LogNormal.(rand(model(x).args.μ, N), rand(model(x).args.σ, N))),
        rand.(LogNormal.(rand(model(x).args.μ, N), rand(model(x).args.σ, N))),
    )'
end

function priorsimulate(x::Induced, N::Int)
    p = priorpredict(x, N)
    rand.(Gamma.(p[1, :], p[2, :]))
end

function posteriorpredict(x::Induced, N::Int)
    n = length(chain(x))
    idxs = sample(1:n, N; replace=true)
    μ₁ = Array(group(chain(x), "μ₁"), [:parameters])[:][idxs]
    μ₂ = Array(group(chain(x), "μ₂"), [:parameters])[:][idxs]
    σ₁ = Array(group(chain(x), "σ₁"), [:parameters])[:][idxs]
    σ₂ = Array(group(chain(x), "σ₂"), [:parameters])[:][idxs]
    hcat(rand.(LogNormal.(μ₁, σ₁)), rand.(LogNormal.(μ₂, σ₂)))'
end

function simulate(x::Induced, N::Int)
    p = posteriorpredict(x, N)
    rand.(Gamma.(p[1, :], p[2, :]))
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

    samples = sample(c.e, 4096)
    @series begin
        seriestype := :histogram
        normalize := true
        xlabel := "y"
        ylabel := "Density"
        label := "Sampled"
        subplot := 4
        samples
    end

    @series begin
        seriestype := :density
        label := "Estimated"
        xlabel := "YFP"
        ylabel := "Density"
        subplot := 4
        linewidth := 3
        xlims := (0, maximum(samples))
        simulate(c, 8192)
    end
end


