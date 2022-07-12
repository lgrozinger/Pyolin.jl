struct PointCharacterisation{T<:Real}
    strain::String63
    backbone::String63
    plasmid::String63
    medians::Vector{T}
end
function PointCharacterisation(strain, backbone, plasmid, df::DataFrame)
    PointCharaterisation(
        String63(strain),
        String63(backbone),
        String63(plasmid),
        median.([Experiment(strain,backbone,plasmid,i,df) for i in IPTGS])
    )
end

const PC = PointCharacterisation

function Base.iterate(x::PC, state=1)
    state > length(x.medians) ? nothing : (x.medians[state], state+1)
end
Base.keys(x::PC)        = IPTGS
Base.length(x::PC)      = length(x.medians)
Base.getindex(x::PC, i) = x.medians[i]

# experimental statistics
StatsBase.mean(x::PC)   = mean.(x)
StatsBase.var(x::PC)    = var.(x)
StatsBase.median(x::PC) = median.(x)


function lsq(X, Y)
    y0, y1 = minimum(Y), maximum(Y)
    x0, x1 = minimum(X), maximum(X)

    @variables x[1:length(X)] k n
    x = collect(x)
    hillexpr = y0 .+ (y1 - y0) .* k^n ./ (x.^n .+ k^n)
    j = Symbolics.jacobian(hillexpr, [k, n])
    J, _ = build_function(j, x, [k, n], expression=Val{false})

    @. F(x, p) = y0 + (y1 - y0) * p[1]^p[2] / (x^p[2] + p[1]^p[2])
    P = [(x1 + x0) / 2, 1.0]
    opts = Dict(
        :upper => [eltype(y)(Inf), eltype(y)(Inf)],
        :lower => [zero(eltype(y)), zero(eltype(y))],
        :showtrace => true
    )
    p = curve_fit(F, J, x, y, P; opts...).params
    y0, y1, p[1], p[2]
end

struct ResponseFunction{T<:Real}
    inputs  :: PointCharacterisation{T}
    outputs :: PointCharacterisation{T}
    ymin :: T
    ymax :: T
    k    :: T
    n    :: T
end
function ResponseFunction(strain, backbone, plasmid, df)
    inputs  = PointCharacterisation(strain, backbone, "1818", df)
    outputs = PointCharacterisation(strain, backbone, plasmid, df)
    y0, y1, k, n = 
end
