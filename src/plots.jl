@recipe function f(e::Experiment)
    (events(e)[CHANNEL], )
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

@userplot CategoryMap
@recipe function f(x::CategoryMap)
    x, y, z = x.args

    categories = reverse(sort(unique(z)))
    function rectangle(x, y)
        [x-1/2, x+1/2, x+1/2, x-1/2, x-1/2], [y-1/2, y-1/2, y+1/2, y+1/2, y-1/2]
    end

    legendposition := :outerright
    xticks := (1:length(x), x)
    yticks := (1:length(y), y)
    xlims := (1/2, length(x) + 1/2)
    ylims := (1/2, length(y) + 1/2)
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
