@recipe function f(e::Experiment)
    (events(e)[CHANNEL], )
end

@recipe function f(x, y::ResponseFunction)
    (x, y.(x))
end

@recipe function f(x, y, z::MultivariateDistribution)
    X = repeat(reshape(x, 1, :), length(y), 1)
    Y = repeat(y, 1, length(x))
    Z = map((x, y) -> pdf(z, [x, y]), X, Y)
    (x, y, Z)
end

@recipe function f(x::AbstractVector{<:Experiment}, y::AbstractVector{<:Experiment})
    (median.(x), median.(y))
end

@userplot StandardisationPlot
@recipe function f(x::StandardisationPlot)
    s, = x.args

    guidefontsize --> 12
    tickfontsize --> 10
    legendfontsize --> 12

    l = @layout [
        [a b c]
        d
    ]
    layout --> l

    @series begin
        seriestype := :density
        xlabel := "k"
        ylabel := "Density"
        label := false
        subplot := 1
        marginalk(s)
    end

    @series begin
        seriestype := :density
        xlabel := L"\theta"
        ylabel := "Density"
        label := false
        subplot := 2
        marginalk(s)
    end

    @series begin
        seriestype := :histogram2d
        xlabel := "k"
        ylabel := L"\theta"
        subplot := 3
        marginalk(s), marginalθ(s)
    end

    @series begin
        seriestype := :histogram
        normalize := true
        xlabel := "y"
        ylabel := "Density"
        label := "Sampled"
        subplot := 4
        sample(s.e, 4096)
    end

    @series begin
        seriestype := :density
        label := "Estimated"
        subplot := 4
        sample(s, 4096)
    end
end

@userplot InverterPlot
@recipe function f(x::InverterPlot)
    if length(x.args) != 2 || !(typeof(x.args[1]) <: AbstractVector{<:Experiment}) ||
        !(typeof(x.args[1]) <: AbstractVector{<:Experiment})
        error("Inverter plots need input and output experiments, not $(typeof(h.args))")
    end

    inputs, outputs = x.args
    fun = ResponseFunction(inputs, outputs)
    xrange = range(minimum(median.(inputs)), maximum(median.(inputs)), 100)

    xlabel --> "Input"
    ylabel --> "Output"
    guidefontsize --> 12
    legendfontsize --> 12
    tickfontsize --> 12
    linewidth --> 3

    @series begin
        seriestype := :line
        label := "Model"
        collect(xrange), fun.(xrange)
    end

    @series begin
        seriestype := :scatter
        markershape := :circle
        markersize := 5
        label := "Data"
        inputs, outputs
    end

    if valid(fun)
        @series begin
            seriestype := :vline
            linewidth := 2
            linestyle := :dash
            label := "Input thresholds"
            [inputlow(fun), inputhigh(fun)]
        end
        @series begin
            seriestype := :hline
            linewidth := 2
            linestyle := :dash
            label := "Output thresholds"
            [outputhigh(fun), outputlow(fun)]
        end
    end
end

@userplot CategoryMap
@recipe function f(x::CategoryMap)
    x, y, z = x.args

    categories = unique(z)
    function rectangle(x, y)
        [x-1/2, x+1/2, x+1/2, x-1/2, x-1/2], [y-1/2, y-1/2, y+1/2, y+1/2, y-1/2]
    end

    legendposition := :outerright
    xticks := (1:length(x), x)
    yticks := (1:length(y), y)
    xlims := (-1/2, length(x) + 1/2)
    ylims := (-1/2, length(y) + 1/2)
    showaxis := :false
    tick_direction := :none
    grid := false

    for category in categories
        @series begin
            seriestype := :shape
            label := string(category)
            xs = []
            ys = []
            for i in 1:length(x)
                for j in 1:length(y)
                    if isequal(category)(z[i, j])
                        a, b = rectangle(i, j)
                        append!(xs, a); push!(xs, NaN)
                        append!(ys, b); push!(ys, NaN)
                    end
                end
            end
            xs, ys
        end
    end
end

@userplot BayesianStandard
@recipe function f(x::BayesianStandard)
    posterior, = x.args

    k = posterior[:, 1]
    θ = posterior[:, 2]
    y = posterior[:, 3]

    Norm = fit(MvNormal, log.(hcat(k, θ))')
    psamples = exp.(rand(Norm, 1024))

    layout := (2,2)

    @series begin
        subplot := 1
        label := false
        seriestype := :scatter
        seriesalpha := 0.25
        xlabel := "k"
        ylabel := "θ"
        psamples[1, :], psamples[2, :]
    end

    @series begin
        subplot := 2
        label := false
        seriestype := :density
        normalize := true
        xlabel := "Expression"
        ylabel := "Density"
        y
    end
    @series begin
        subplot := 3
        label := false
        seriestype := :density
        normalize := true
        xlabel := "k (mRNA per cycle)"
        ylabel := "Density"
        k
    end
    @series begin
        subplot := 4
        label := false
        seriestype := :density
        normalize := true
        xlabel := "θ (protein per mRNA)"
        ylabel := "Density"
        θ
    end
end
