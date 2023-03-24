function transfer(x, p)
    y₀, y₁, k, n = p
    return y₀ .+ y₁ ./ (1 .+ (x ./ k).^n)
end

function lsq(X::Vector{T}, Y::Vector{T}, W::Vector{T}) where {T<:Real}
    y₀, y₁ = minimum(Y), maximum(Y)
    x₀, x₁ = minimum(X), maximum(X)
    opts = Dict(
        :upper => [T(Inf), T(Inf), T(Inf), T(Inf)],
        :lower => [eps(T), eps(T), eps(T), one(T)],
        :show_trace => false,
        :maxIter => 1024,
        :autodiff => :forwarddiff,
    )
    return curve_fit(transfer, X, Y, W, [y₀, y₁, x₀/2 + x₁/2, one(T)]; opts...)
end

function lsq(X::Vector{T}, Y::Vector{T}) where {T<:Real}
    y₀, y₁ = minimum(Y), maximum(Y)
    x₀, x₁ = minimum(X), maximum(X)
    
    opts = Dict(
        :upper => [T(Inf), T(Inf), T(Inf), T(Inf)],
        :lower => [eps(T), eps(T), eps(T), eps(T)],
        :show_trace => false,
        :maxIter => 1024,
        :autodiff => :forwarddiff,
    )
    return curve_fit(transfer, X, Y, [y₀, y₁, x₀/2 + x₁/2, one(T)]; opts...)
end

struct Hill{T<:Real}
    strain::String
    backbone::String
    plasmid::String
    fit::LsqFit.LsqFitResult
    ins::Vector{T}
    outs::Vector{T}
end

function Hill(inputs::Vector{<:Experiment}, outputs::Vector{<:Experiment})
    s, b, p = outputs[1].strain, outputs[1].backbone, outputs[1].plasmid
    function nonzero(x::Vector{T}) where {T<:Real}
        min = minimum(x)
        if min < 0
            return x .- min
        else
            return x
        end
    end
    X, Y = nonzero(median.(inputs)), nonzero(median.(outputs))
    W = 1 ./ var.(outputs)
    fit = lsq(X, Y, W)
    Hill(s, b, p, fit, X, Y)
end

function Hill(strain, backbone, plasmid)
    inputs  = Experiments(search(strain, backbone, INPUT))
    autos   = Experiments(search(strain, backbone, AUTOFLUOR))
    stds    = Experiments(search(strain, backbone, STANDARD))
    outputs = Experiments(search(strain, backbone, plasmid))
    Hill(rpuconvert(inputs, autos, stds), rpuconvert(outputs, autos, stds))
end

Hills(fn::AbstractString) = Hills(CSV.read(fn, DataFrame))

y0(R::Hill) = R.fit.param[1]
y1(R::Hill) = R.fit.param[2] + y0(R)
K(R::Hill)  = R.fit.param[3]
N(R::Hill)  = R.fit.param[4]

function (R::Hill)(x::T) where {T<:Real}
    y0(R) + (y1(R) - y0(R)) * K(R)^N(R) / (x^N(R) + K(R)^N(R))
end

function inputs(x::Hill)
    Experiments(search(x.strain, x.backbone, INPUT))
end
function outputs(x::Hill)
    Experiments(search(x.strain, x.backbone, x.plasmid))
end
function autos(x::Hill)
    Experiments(search(x.strain, x.backbone, AUTOFLUOR))
end
function standards(x::Hill)
    Experiments(search(x.strain, x.backbone, STANDARD))
end

repressor(x::Hill)              = split(x.plasmid, "_")[1]
samestrain(x::Hill, y::Hill)    = x.strain == y.strain
samerepressor(x::Hill, y::Hill) = repressor(x) == repressor(y)
orthogonal(x::Hill, y::Hill)    = !samestrain(x, y) || !samerepressor(x, y)
outputhigh(x::Hill)             = x.ymax / 2
outputlow(x::Hill)              = x.ymin * 2
function inputlow(x::Hill)
    (K(x)^N(x) * (2 / (1 - y0(x) / (y1(x) - y0(x))) - 1))^(1/N(x))
end
function inputhigh(x::Hill)
    y₀, y₁, k, n = x.ymin, x.ymax, x.K, x.n
    (K(x)^N(x) * (y1(x) / y0(x) - 2))^(1/N(x))
end
valid(x::Hill)                  = outputhigh(x) > outputlow(x) #&& inputhigh(x) > inputlow(x)
invertshigh(x::Hill, y::Hill)   = y(outputhigh(x)) < outputlow(y)
invertslow(x::Hill, y::Hill)    = y(outputlow(x)) > outputhigh(y)
inverts(x::Hill, y::Hill)       = invertshigh(x, y) && invertslow(x, y)
compatible(x::Hill, y::Hill)    = orthogonal(x, y) && valid(x) && valid(y) && inverts(x, y)

istable(x::Vector{T}) where {T<:Hill}   = true
rowaccess(x::Vector{T}) where {T<:Hill} = true
rows(x::Vector{T}) where {T<:Hill}      = x
schema(x::Vector{T}) where {T<:Hill}    = Tables.Schema(fieldnames(T), fieldtypes(T))

@recipe function f(x, y::Hill)
    (x, y.(x))
end

@userplot NotPlot
@recipe function f(x::NotPlot)
    h, = x.args
    xs = range(0.8*minimum(h.ins), 1.125*maximum(h.ins), 100)

    xlabel --> "Input"
    ylabel --> "Output"
    guidefontsize  --> 10
    legendfontsize --> 10
    tickfontsize   --> 10
    palette --> :Accent_8

    @series begin
        seriestype := :line
        linewidth --> 3
        label := "Hill model"
        collect(xs), h.(xs)
    end

    @series begin
        seriestype := :scatter
        markershape := :circle
        markersize --> 5
        linewidth := 1
        label := "Data points"
        h.ins, h.outs
    end

    if valid(h) && get(plotattributes, :thresholds, false)
        @series begin
            seriestype := :vline
            linewidth := 2
            linestyle := :dash
            label := "Input thresholds"
            [inputlow(h), inputhigh(h)]
        end
        @series begin
            seriestype := :hline
            linewidth := 2
            linestyle := :dash
            label := "Output thresholds"
            [outputhigh(h), outputlow(h)]
        end
    end
end

@userplot NotNotPlot
@recipe function f(x::NotNotPlot)
    h, = x.args
    xs = range(0.8*minimum(h.ins), 1.125*maximum(h.ins), 100)

    xlabel --> "Output"
    ylabel --> "Input"
    guidefontsize  --> 10
    legendfontsize --> 10
    tickfontsize   --> 10
    palette --> :Accent_8

    @series begin
        seriestype := :line
        linewidth --> 3
        label := "Hill model"
        h.(xs), collect(xs)
    end

    @series begin
        seriestype := :scatter
        markershape := :circle
        markersize --> 5
        linewidth := 1
        label := "Data points"
        h.outs, h.ins
    end

    if valid(h) && get(plotattributes, :thresholds, false)
        @series begin
            seriestype := :hline
            linewidth := 2
            linestyle := :dash
            label := "Input thresholds"
            [inputlow(h), inputhigh(h)]
        end
        @series begin
            seriestype := :vline
            linewidth := 2
            linestyle := :dash
            label := "Output thresholds"
            [outputhigh(h), outputlow(h)]
        end
    end
end
