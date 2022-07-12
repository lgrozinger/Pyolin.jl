"""A gamma distribution mixed model"""
@model function model1717(y::Vector{Float64}, μ::Vector{Float64}, σ::Matrix{Float64})
    p ~ MvLogNormal(μ, σ)
    y ~ filldist(Gamma(p[1], p[2]), length(y))
    return y
end
@model function model1717(N::Int, μ::Vector{Float64}, σ::Matrix{Float64})
    y = Vector{Float64}(undef, N)
    p ~ MvLogNormal(μ, σ)
    y ~ filldist(Gamma(p[1], p[2]), N)
    return y
end

function model1717(e::Experiment, N::Int, μ::Vector{Float64}, σ::Matrix{Float64}; witherror=true)
    y = Float64.(sample(events(e)[CHANNEL], N; replace=false))
    if witherror
        a = Experiment(e.strain, e.backbone, "1201", e.iptg)
        model1717(y, μ, σ, fit(Normal{Float64}, a))
    else
        model1717(y, μ, σ)
    end
end

"""A gamma distribution mixed model with dose-reponse"""
@model function model1818(x, y, μ::Vector{T}, σ::Matrix{T}) where {T<:Real}
    if y === missing
        y = Vector{T}(undef, length(x))
    end
    p ~ MvLogNormal(μ, σ)
    k = p[1] .* x.^p[4] ./ (x.^p[4] .+ p[3]^p[4])
    y .~ Gamma.(k, p[2])
    return y
end

function cond1818(es::Vector, N::Int, μ, σ)
    x = reduce(vcat, repeat([i], N) for i in getproperty.(es, :iptg))
    y = reduce(vcat, sample(events(e)[CHANNEL], N; replace=false) for e in es)
    model1818(x, y, μ, σ)
end

function StandardModel(model, N::Int)
    sampler = NUTS(256, 0.65)
    Nt = Threads.nthreads()
    opts = Dict(
        :discard_adapt => true,
        :progress => false
    )
    chain = sample(model, sampler, MCMCThreads(), N, Nt; opts...)
    psamples = reduce(hcat, x[:] for x in get_params(chain).p)
    mvn = fit(MvNormal, log.(psamples)')
    pdist = MvLogNormal(mvn.μ, collect(mvn.Σ))
    N -> reduce(hcat, posterior() for _ in 1:N)
end

function InputModel(model, N::Int)
    sampler = NUTS(256, 0.65)
    Nt = Threads.nthreads()
    chain = sample(model, sampler, MCMCThreads(), N, Nt)
    psamples = reduce(hcat, x[:] for x in get_params(chain).p)
    mvn = fit(MvNormal, log.(psamples)')
    model1818(model.args.x, missing, mvn.μ, collect(mvn.Σ))
end
