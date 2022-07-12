struct ResponseFunction{T<:Real}
    strain::String63
    backbone::String63
    plasmid::String63
    ymin::T
    ymax::T
    K::T
    n::T
end
function ResponseFunction(inputs::Vector{<:Experiment}, outputs::Vector{<:Experiment})
    strain   = outputs[1].strain
    backbone = outputs[1].backbone
    plasmid  = outputs[1].plasmid

    X, Y = median.(inputs), median.(outputs)
    ymin, ymax, K, n = lsqfitting(X, Y)
    ResponseFunction(strain, backbone, plasmid, ymin, ymax, K, n)
end
function ResponseFunction(strain, backbone, plasmid)
    inputs  = Experiments(search(strain, backbone, INPUT))
    outputs = Experiments(search(strain, backbone, plasmid))
    ResponseFunction(inputs, outputs)
end
function ResponseFunctions(fn::AbstractString)
    function rowfun(r)
        ResponseFunction(
            String63(r.strain), String63(r.backbone), String63(r.plasmid),
            r.ymin, r.ymax, r.K, r.n
        )
    end
    rowfun.(eachrow(CSV.read(fn, DataFrame)))
end

const RF = ResponseFunction

# the function 
function (R::RF)(x::T) where {T<:Real}
    y0, y1, n, K = R.ymin, R.ymax, R.n, R.K
    y0 + (y1 - y0) * K^n / (x^n + K^n)
end

# properties of response functions
repressor(x::RF) = split(x.plasmid, "_")[1]
samestrain(x::RF, y::RF) = x.strain == y.strain
samerepressor(x::RF, y::RF) = repressor(x) == repressor(y)
orthogonal(x::RF, y::RF) = !samestrain(x, y) || !samerepressor(x, y)
outputhigh(x::RF) = x.ymax / 2
outputlow(x::RF) = x.ymin * 2
inputhigh(x::RF) = (x.K^x.n * (x.ymax - 2 * x.ymin) / x.ymin)^(1 / x.n)
inputlow(x::RF)  = (x.K^x.n * x.ymax / (x.ymax - 2 * x.ymin))^(1 / x.n)
valid(x::RF) = outputhigh(x) > outputlow(x) && inputhigh(x) > inputlow(x)
invertshigh(x::RF, y::RF) = y(outputhigh(x)) < outputlow(y)
invertslow(x::RF, y::RF) = y(outputlow(x)) > outputhigh(y)
inverts(x::RF, y::RF) = invertshigh(x, y) && invertslow(x, y)
compatible(x::RF, y::RF) = orthogonal(x, y) && valid(x) && valid(y) && inverts(x, y)

# tables interface
istable(x::Vector{T}) where {T<:RF}   = true
rowaccess(x::Vector{T}) where {T<:RF} = true
rows(x::Vector{T}) where {T<:RF}      = x
schema(x::Vector{T}) where {T<:RF}    = Tables.Schema(fieldnames(T), fieldtypes(T))


