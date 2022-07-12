using Pyolin
using HypothesisTests
using StatsBase

function KSTest(E::Experiment, F::Experiment, n::Int)
    e = sample(events(E), n; replace=false)
    f = sample(events(F), n; replace=false)
    ApproximateTwoSampleKSTest(e, f)
end

function KSTest(E::Experiment, F::Experiment)
    n, m = length(events(E)), length(events(F))
    n < m ? KSTest(E, F, n) : KSTest(E, F, m)
end

function KSTest(A::Gate, B::Gate)
    n = length(A.outputs)
    results = Vector{ApproximateTwoSampleKSTest}(undef, n)
    for i in eachindex(collect(zip(A.outputs, B.outputs)))
        results[i] = KSTest(A.outputs[i], B.outputs[i])
    end
    results
end

function TVTest(E::Experiment, F::Experiment, n::Int)
    e = sample(events(E), n; replace=false)
    f = sample(events(F), n; replace=false)
    totalvariation(e, f)
end

function TVTest(E::Experiment, F::Experiment)
    n, m = length(events(E)), length(events(F))
    n < m ? TVTest(E, F, n) : TVTest(E, F, m)
end

function TVTest(A::Gate, B::Gate)
    n = length(A.outputs)
    results = Vector{Float32}(undef, n)
    for i in eachindex(collect(zip(A.outputs, B.outputs)))
        results[i] = TVTest(A.outputs[i], B.outputs[i])
    end
    results
end

a = Gate(strain="CC118Lpir", backbone="pSeva221", plasmid="1201", rpu=false)
b = Gate(strain="CC118Lpir", backbone="pSeva231", plasmid="1201", rpu=false)


