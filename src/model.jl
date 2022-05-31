using Turing
using ReverseDiff
using MCMCChains
using FillArrays

N(μ, σ₁, σ₂, N) = rand.(Normal.(rand(Normal(μ, σ₁), N), σ₂))
R(μ, σ, N) = rand.(Poisson.(exp.(rand(Normal(μ, σ), N))))

@model function transcription(r)
    n = length(r)
    μ₁ ~ Normal()
    σ₁ ~ truncated(Normal(), 0, Inf)

    λ ~ filldist(Normal(μ₁, σ₁), length(r))
    for i = 1:length(r)
        r[i] ~ Poisson(exp(λ[i]))
    end
end

model = transcription(R(0, 1, 2^10))
sampler = NUTS(256, 0.9)
chain = sample(model, sampler, 2^12; discard_adapt=true)

# samples = N(2, 1, 2, 2^18)
# posterior_samples = N(mean(chain[:μ]), var(chain[:μ]), mean(chain[:σ]), 2^18)



# setprogress!(false)
# model = m(r)
# sampler = NUTS(2^8, 0.65)
# nchains = Threads.nthreads()
# chain = sample(model, sampler, MCMCThreads(), 2^10, nchains; discard_adapt=true)
