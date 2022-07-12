const DEFGATE! = (Spike!("FSC-H", 30, 4, 0.1) ∘ NonZero!("FSC-A", "SSC-A", "B1-H", "B1-A"))

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
Experiments(frame::DataFrame) = Experiment.(eachrow(frame))
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

# experimental data
row(e::Experiment) = search(e.strain, e.backbone, e.plasmid, e.iptg)
events(e::Experiment) = DEFGATE!(events(row(e)))
StatsBase.median(e::Experiment) = e.median
StatsBase.mean(e::Experiment) = e.mean
StatsBase.var(e::Experiment) = e.var
StatsBase.maximum(e::Experiment) = e.max
StatsBase.minimum(e::Experiment) = e.min

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
function rpuconvert(x::Experiment{T}, auto::Experiment{T}, standard::Experiment{T}) where {T<:Real}
    if median(x) < median(auto)
        @warn "RPU conversion smudge for $(x.strain), $(x.backbone), $(x.plasmid), $(x.iptg)\nMedian of experiment: $(x.median)\nMedian of autofluorescence: $(auto.median)"
        c = eps(T)^2
    else
        c = (median(x) - median(auto)) / median(x) / (median(standard) - median(auto))
    end
    Experiment(
        x.strain, x.backbone, x.plasmid, x.iptg,
        x.mean * c, x.median * c, x.var * c^2, x.max * c, x.min * c,
        x.fn
    )
end
rpuconvert(x::T, auto::T, standard::T) where {T<:Vector{<:Experiment}} = rpuconvert.(x, auto, standard)
