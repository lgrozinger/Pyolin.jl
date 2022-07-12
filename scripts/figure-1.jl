using Pyolin
using Distributions
using StatsPlots
using CSV
using DataFrames
using Plots.PlotMeasures


include("figure_settings.jl")

function figure1_1(;strain, backbone, iptg, fn, kwargs...)
    auto = Experiment(strain=strain, backbone=backbone, plasmid="1201", iptg=iptg)
    standard = Experiment(strain=strain, backbone=backbone, plasmid="1717", iptg=iptg)
    plt = histogram(
        events(auto);
        DEFS...,
        label="Autofluorescence",
        color=RGBA(153/255, 153/255, 153/255, 0.53),
        ylabel="Sample frequency",
        xlabel="Log intensity (channel B1-H A.U.)",
        size = (450, 300),
        kwargs...
            )

    histogram!(
        plt,
        events(standard);
        color=RGBA(1, 1, 0, 0.53),
        DEFS...,
        label="Standardisation",
        kwargs...
            )

    dist = fit(Gamma, events(standard))
    x0 = minimum([minimum(auto), minimum(standard)])
    mx = x0
    while cdf(dist, mx) < 0.9995
        mx = mx + 0.01f0
    end

    plot!(plt; xlims=(x0, mx), xscale=:log10, legend=:best, bottom_margin=5mm, kwargs...)
    savefig(plt, fn)
    plt
end

#figure 1.1
figure1_1(
    ;strain="KT2440",
    backbone="pSeva231",
    iptg=0,
    fn="/home/lewis/thesis/figures/figure1.1.svg"
);




