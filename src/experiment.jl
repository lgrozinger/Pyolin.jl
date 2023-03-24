using FileIO
using StatsBase
using Distributions

const DEFGATE! = (Zscore!("B1-A", 3.5) ∘ Zscore!("SSC-H", 3.5) ∘ Zscore!("FSC-H", 3.5) ∘ NonZero!("FSC-A", "SSC-A", "B1-H", "B1-A"))

abstract type Experiment{T<:Real} end

function (::Type{T})(x::DataFrameRow) where {T<:Experiment}
    flowsample = FileIO.load(DATADIR * x.filename)
    return T(
        convert(String, x.strain),
        convert(String, x.backbone),
        convert(String, x.plasmid),
        x.iptg,
        flowsample
    )
end

function (::Type{T})(x::DataFrame) where {T<:Experiment}
    return T.(eachrow(x))
end

function for_other_plasmid(e::Experiment, other_plasmid)
    match(x) = all(getproperty(x, s) == getproperty(e, s) for s in [:strain, :backbone, :iptg])
    return filter(x -> x.plasmid == other_plasmid, filter(match, index))
end

Auto(e::T) where {T<:Experiment}     = T(for_other_plasmid(e, "1201"))[1]
Input(e::T) where {T<:Experiment}    = T(for_other_plasmid(e, "1818"))[1]
Standard(e::T) where {T<:Experiment} = T(for_other_plasmid(e, "1717"))[1]

Autos(E)     = Auto.(E)
Inputs(E)    = Input.(E)
Standards(E) = Standard.(E)

Base.size(e::Experiment, args...) = size(events(e), args...)
StatsBase.median(e::Experiment)  = median(events(e))
StatsBase.mean(e::Experiment)    = mean(events(e))
StatsBase.var(e::Experiment)     = var(events(e))
StatsBase.maximum(e::Experiment) = maximum(events(e))
StatsBase.minimum(e::Experiment) = minimum(events(e))

function StatsBase.sample(e::Experiment, N::Int; replace=false)
    return sample(events(e), N; replace=replace)
end

function Distributions.fit(d::Type{Normal}, e::Experiment)
    return Normal(mean(e), sqrt(var(e)))
end

function Distributions.fit(d::Type{LogNormal}, e::Experiment)
    return LogNormal(mean(e), sqrt(var(e)))
end

function Distributions.fit(d::Type{<:UnivariateDistribution}, e::Experiment)
    return fit(d, events(e)[CHANNEL])
end

struct RawExperiment{T<:Real, I} <: Experiment{T}
    strain::String
    backbone::String
    plasmid::String
    iptg::Int
    flow::FlowSample{T, I}
end

function events(e::RawExperiment)
    return Vector(DEFGATE!(e.flow).data[CHANNEL])
end

struct RPUExperiment{T<:Real, I} <: Experiment{T}
    strain::String
    backbone::String
    plasmid::String
    iptg::Int
    flow::FlowSample{T, I}
end

function events(e::RPUExperiment)
    raw = RawExperiment(e.strain, e.backbone, e.plasmid, e.iptg, e.flow)
    auto = median(Auto(raw))
    standard = median(Standard(raw))
    return events(raw) ./ standard
#    return (events(raw) .- auto) ./ (standard - auto)
end
