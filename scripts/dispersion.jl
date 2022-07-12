using Pyolin
using DataFrames
using ProgressMeter
using CSV
using LaTeXStrings
using Plots
using Plots.PlotMeasures

include("figure_settings.jl")
FIGDIR = "/home/lewis/thesis/hostchapter/figures/dispersion/"

statsframe = CSV.read("/home/lewis/sauce/julia/Pyolin/data/experimentstats.csv", DataFrame)

function clustered_scatter(df, groups, fx, fy; kwargs...)
    plt = plot(;kwargs...)
    dfs = groupby(df, groups)
    for group in dfs
        label = reduce(*, [""; [string(getproperty(group, k)[1]) * " " for k in groups]])
        scatter!(
            plt,
            fx.(eachrow(group)),
            fy.(eachrow(group));
            label = label[1:end-1],
            kwargs...
            )
    end
    plt
end

autos = filter(:plasmid => ==("1201"), statsframe)
standards = filter(:plasmid => ==("1717"), statsframe)
inputs = filter(:plasmid => ==("1818"), statsframe)
gates = filter(:plasmid => ∉(["1201", "1818", "1717"]), statsframe)

automeans = clustered_scatter(
    autos,
    [:strain, :backbone],
    r -> r.shape,
    r -> r.mean;
    legend = :bottomright,
    ylabel = "Sample mean",
    xlabel = "Model mean",
)
savefig(automeans, FIGDIR * "automeans.svg")

autovars = clustered_scatter(
    autos,
    [:strain, :backbone],
    r -> r.scale^2,
    r -> r.var;
    legend = :bottomright,
    ylabel = "Sample variance",
    xlabel = "Model variance",
    xscale = :log10,
    yscale = :log10
)
savefig(autovars, FIGDIR * "autovars.svg")

stdmeans = clustered_scatter(
    standards,
    [:strain, :backbone],
    r -> r.shape * r.scale,
    r -> r.mean;
    xscale = :log10,
    yscale = :log10,
    legend = :bottomright,
    ylabel = "Sample mean",
    xlabel = "Model mean"
)
plot!(stdmeans, [1e0, 1e2], [1e0, 1e2], label="")
savefig(stdmeans, FIGDIR * "stdmeans.svg")

stdvars = clustered_scatter(
    standards,
    [:strain, :backbone],
    r -> r.shape * r.scale ^ 2,
    r -> r.var;
    xscale = :log10,
    yscale = :log10,
    legend = :bottomright,
    ylabel = "Sample variance",
    xlabel = "Model variance"
)
plot!(stdvars, [1e-1, 1e3], [1e-1, 1e3], label="")
savefig(stdvars, FIGDIR * "stdvars.svg")

inputmeans = clustered_scatter(
    inputs,
    [:strain, :backbone],
    r -> r.shape * r.scale,
    r -> r.mean;
    xscale = :log10,
    yscale = :log10,
    legend = :bottomright,
    ylabel = "Sample mean",
    xlabel = "Model mean"
)
plot!(inputmeans, [5e-1, 5e1], [5e-1, 5e1], label="")
savefig(inputmeans, FIGDIR * "inputmeans.svg")

inputvars = clustered_scatter(
    inputs,
    [:iptg],
    r -> r.shape * r.scale ^ 2,
    r -> r.var;
    xscale = :log10,
    yscale = :log10,
    legend = :bottomright,
    ylabel = "Sample variance",
    xlabel = "Model variance"
)
plot!(inputvars, [5e-2, 1e2], [5e-2, 1e2], label="")
savefig(inputvars, FIGDIR * "inputvars.svg")

gatemeans = clustered_scatter(
    gates,
    [],
    r -> r.shape * r.scale,
    r -> r.mean;
    xscale = :log10,
    yscale = :log10,
    legend = :bottomright,
    ylabel = "Sample mean",
    xlabel = "Model mean"
)
plot!(gatemeans, [5e-1, 5e2], [5e-1, 5e2], label="")
savefig(gatemeans, FIGDIR * "gatemeans.svg")

gatevars = clustered_scatter(
    gates,
    [],
    r -> r.shape * r.scale ^ 2,
    r -> r.var;
    xscale = :log10,
    yscale = :log10,
    legend = :bottomright,
    ylabel = "Sample variance",
    xlabel = "Model variance"
)
plot!(gatevars, [1e-2, 1e4], [1e-2, 1e4], label="")
savefig(gatevars, FIGDIR * "gatevars.svg")


function bayesianstats_standard(row)
    e = Experiment(row.strain, row.backbone, row.plasmid, row.iptg)
    m = proposemodel(e, 4096, [1.0, 1.0] , [0.1 0.0; 0.0 0.1])
    s = reduce(vcat, m() for _ in 1:2^15)
    mean(s), var(s)
end

standard_stats = [tuple(r.mean, r.var, bayesianstats_standard(r)...) for r in eachrow(standards)]

meanplt = scatter(
    getindex.(standard_stats, 1),
    getindex.(standard_stats, 3);
    xlabel = "Sample mean",
    ylabel = "Model mean",
    xscale = :log10,
    yscale = :log10
)
plot!(meanplt, [1e0, 1e2], [1e0, 1e2], label=false)

varplt = scatter(
    getindex.(standard_stats, 4),
    getindex.(standard_stats, 2);
    xlabel = "Model variance",
    ylabel = "Sample variance",
    xscale = :log10,
    yscale = :log10
)
plot!(varplt, [1e-1, 1e3], [1e-1, 1e3], label=false)
