const DEFGATE! = (Zscore!("B1-H", 3.5) ∘ Zscore!("SSC-H", 3.5) ∘ Zscore!("FSC-H", 3.5) ∘ NonZero!("FSC-A", "SSC-A", "B1-H", "B1-A"))

struct Experiment{T<:Real}
    strain::String63
    backbone::String63
    plasmid::String63
    iptg::Int
    mean::T
    median::T
    var::T
    max::T
    min::T
    fn::String255
end

function Experiment(strain, backbone, plasmid, iptg)
    row = search(strain, backbone, plasmid, iptg)
    fn = String255(row.filename)
    samples = DEFGATE!(events(row))

    Experiment(
        String63(strain),
        String63(backbone),
        String63(plasmid),
        iptg,
        mean(samples[CHANNEL]),
        median(samples[CHANNEL]),
        var(samples[CHANNEL]),
        maximum(samples[CHANNEL]),
        minimum(samples[CHANNEL]),
        fn
    )
end

Experiment(row::DataFrameRow) = Experiment(row.strain, row.backbone, row.plasmid, row.iptg)

Experiments(frame) = Experiment.(eachrow(frame))

function Experiments(fn::AbstractString)
    function rowfun(r)
        Experiment(
            String63(r.strain), String63(r.backbone), String63(r.plasmid), r.iptg,
            r.mean, r.median, r.var, r.max, r.min, String255(r.fn)
        )
    end
    rowfun.(eachrow(CSV.read(fn, DataFrame)))
end

function Experiments(strain, backbone, plasmid, fn::AbstractString)
    function rowfun(r)
        Experiment(
            String63(r.strain), String63(r.backbone), String63(r.plasmid), r.iptg,
            r.mean, r.median, r.var, r.max, r.min, String255(r.fn)
        )
    end
    frame = CSV.read(fn, DataFrame)
    filter!(r -> r.strain == strain && r.backbone == backbone && r.plasmid == plasmid, frame)
    sort!(frame, [:iptg])
    rowfun.(eachrow(frame))
end

function Experiments(strain, backbone, plasmid)
    iptgs = sort(unique(index.iptg))
    [Experiment(strain, backbone, plasmid, iptg) for iptg in iptgs]
end

Auto(e::Experiment)     = Experiment(e.strain, e.backbone, String63("1201"), e.iptg)
Standard(e::Experiment) = Experiment(e.strain, e.backbone, String63("1717"), e.iptg)
Input(e::Experiment)    = Experiment(e.strain, e.backbone, String63("1818"), e.iptg)
Autos(E::Vector{<:Experiment})     = Auto.(E)
Standards(E::Vector{<:Experiment}) = Standard.(E)
Inputs(E::Vector{<:Experiment}) = Input.(E)


# experimental data
row(e::Experiment) = search(e.strain, e.backbone, e.plasmid, e.iptg)
events(e::Experiment) = DEFGATE!(events(row(e)))
Base.length(e::Experiment) = length(events(e)[CHANNEL])
StatsBase.median(e::Experiment) = e.median
StatsBase.mean(e::Experiment) = e.mean
StatsBase.var(e::Experiment) = e.var
StatsBase.maximum(e::Experiment) = e.max
StatsBase.minimum(e::Experiment) = e.min
StatsBase.quantile(e::Experiment, α::T) where {T<:Real} = quantile(events(e)[CHANNEL], α)
StatsBase.sample(e::Experiment, N::Int; replace=false) = sample(events(e)[CHANNEL], N; replace=replace)

# fitting
Distributions.fit(d::Type{Normal}, e::Experiment) = d(e.mean, √(e.var))
function Distributions.fit(d::Type{<:UnivariateDistribution}, e::Experiment; g=DEFGATE!)
    d(params(fit(d, g(events(e))[CHANNEL]))...)
end

# tables interface
columnnames(e::Experiment) = schema([e]).names
getcolumn(e::Experiment, i::Int) = getproperty(e, columnnames(e)[i])
getcolumn(e::Experiment, s::Symbol) = getproperty(e, s)

istable(x::Vector{Experiment})             = true
rowaccess(x::Vector{Experiment})           = true
rows(x::Vector{Experiment})                = x
schema(x::Vector{T}) where {T<:Experiment} = Tables.Schema(fieldnames(T), fieldtypes(T))

# rpu conversion
function rpuevents(x::Experiment{T}, standard::Experiment{T}) where {T<:Real}
    c = 1 / median(standard)
    sample(x, length(x)) .* c
end
function rpuevents(x::Experiment{T}, auto::Experiment{T}, standard::Experiment{T}) where {T<:Real}
    if median(x) < median(auto)
        @warn "RPU conversion smudge for $(x.strain), $(x.backbone), $(x.plasmid), $(x.iptg)\nMedian of experiment: $(x.median)\nMedian of autofluorescence: $(auto.median)"
        c = eps(T)
    else
        c = (median(x) - median(auto)) / median(x) / (median(standard) - median(auto))
    end
    sample(x, length(x)) .* c
end
rpuevents(x::AbstractVector, auto::AbstractVector, standard::AbstractVector) = rpuevents.(x, auto, standard)
rpuevents(x::AbstractVector, standard::AbstractVector) = rpuevents.(x, standard)

function rpuconvert(x::Experiment, auto::Experiment, standard::Experiment)
    if median(x) < median(auto)
        @warn "RPU conversion smudge for $(x.strain), $(x.backbone), $(x.plasmid), $(x.iptg)\nMedian of experiment: $(x.median)\nMedian of autofluorescence: $(auto.median)"
        c = eps(typeof(x).parameters[1])
    else
        c = (median(x) - median(auto)) / median(x) / (median(standard) - median(auto))
    end
    Experiment(
        x.strain, x.backbone, x.plasmid, x.iptg,
        x.mean * c, x.median * c, x.var * c^2, x.max * c, x.min * c,
        x.fn
    )
end
rpuconvert(x::AbstractVector, auto::AbstractVector, standard::AbstractVector) = rpuconvert.(x, auto, standard)

function rpuconvert(x::Experiment, standard::Experiment)
    c = 1 / median(standard)
    Experiment(
        x.strain, x.backbone, x.plasmid, x.iptg,
        x.mean * c, x.median * c, x.var * c^2, x.max * c, x.min * c,
        x.fn
    )
end
rpuconvert(x::AbstractVector, standard::AbstractVector) = rpuconvert.(x, standard)

