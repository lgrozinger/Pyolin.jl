function _f1_(y₀, y₁)
    function (x, p)
        k, n = p
        ((y₁ - y₀) * k^n) ./ (k^n .+ x.^n) .+ y₀
    end
end
function _j1_(y₀, y₁, N)
    @variables k n x[1:N]
    x = collect(x)
    expression = ((y₁ - y₀) * k^n) ./ (k^n .+ x.^n) .+ y₀
    j = Symbolics.jacobian(expression, [k, n])
    J, _ = build_function(j, x, [k, n], expression=Val{false})
    J
end

function lsq(X::Vector{T}, Y::Vector{T}) where {T<:Real}
    y₀, y₁ = minimum(Y), maximum(Y)
    x₀, x₁ = minimum(X), maximum(X)
    opts = Dict(
        :upper => [T(Inf), T(256.0)],
        :lower => [x₀, one(T)],
        :show_trace => false,
        :maxIter => 1000
    )
    try
        k, n = curve_fit(
            _f1_(y₀, y₁),
            _j1_(y₀, y₁, length(X)),
            X,
            Y,
            [x₀/2 + x₁/2, one(T)];
            opts...,
        ).param
        y₀, y₁, k, n
    catch e
        if e isa LinearAlgebra.SingularException
            T(0.0), T(0.0), T(0.0), T(1.0)
        else
            throw(e)
        end
    end
end

struct Hill{T<:Real}
    strain::String63
    backbone::String63
    plasmid::String63
    ymin::T
    ymax::T
    K::T
    n::T
    ins::Vector{T}
    outs::Vector{T}
end

function Hill(inputs::Vector{<:Experiment}, outputs::Vector{<:Experiment})
    s, b, p = outputs[1].strain, outputs[1].backbone, outputs[1].plasmid
    X, Y = median.(inputs), median.(outputs)
    y₀, y₁, k, n = lsq(X, Y)
    Hill(s, b, p, y₀, y₁, k, n, X, Y)
end

function Hill(strain, backbone, plasmid)
    inputs  = Experiments(search(strain, backbone, INPUT))
    autos   = Experiments(search(strain, backbone, AUTOFLUOR))
    stds    = Experiments(search(strain, backbone, STANDARD))
    outputs = Experiments(search(strain, backbone, plasmid))
    Hill(rpuconvert(inputs, autos, stds), rpuconvert(outputs, autos, stds))
end

function Hills(df)
    function rowfun(r)
        Hill(
            String63(r.strain), String63(r.backbone), String63(r.plasmid),
            r.ymin, r.ymax, r.K, r.n
        )
    end
    rowfun.(eachrow(df))
end

Hills(fn::AbstractString) = Hills(CSV.read(fn, DataFrame))

function (R::Hill)(x::T) where {T<:Real}
    y0, y1, n, K = R.ymin, R.ymax, R.n, R.K
    y0 + (y1 - y0) * K^n / (x^n + K^n)
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
    y₀, y₁, k, n = x.ymin, x.ymax, x.K, x.n
    (k^n * (2 / (1 - y₀ / (y₁ - y₀)) - 1))^(1/n)
end
function inputhigh(x::Hill)
    y₀, y₁, k, n = x.ymin, x.ymax, x.K, x.n
    (k^n * (y₁ / y₀ - 2))^(1/n)
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
