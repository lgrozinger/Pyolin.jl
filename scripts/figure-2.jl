using Pyolin
using Distributions
using StatsPlots
using Plots.PlotMeasures


include("figure_settings.jl")

function figure1_2(;strain, backbone, iptg, fn, kwargs...)
    auto = Experiment(strain=strain, backbone=backbone, plasmid="1201", iptg=iptg)
    standard = Experiment(strain=strain, backbone=backbone, plasmid="1717", iptg=iptg)
    inrpu= Experiment(
        strain=strain,
        backbone=backbone,
        plasmid="1818",
        iptg=iptg,
        auto=auto,
        standard=standard
    )
    inau= Experiment(
        strain=strain,
        backbone=backbone,
        plasmid="1818",
        iptg=iptg
    )

    audist = fit(Gamma, events(inau))
    rpudist = fit(Gamma, events(inrpu))
    x0 = minimum([events(inau); events(inrpu)])
    mx = x0
    while cdf(audist, mx) < 0.9995 || cdf(rpudist, mx) < 0.9995
        mx = mx + 0.01f0
    end

    @show x0, mx

    plt = histogram(
        filter(x -> x0 < x < mx, events(inrpu));
        DEFS...,
        label="RPU",
        ylabel="Sample frequency",
        xlabel="Log intensity (channel B1-H A.U.)",
        size = (600, 350),
        kwargs...
            )

    histogram!(
        plt,
        filter(x -> x0 < x < mx, events(inau));
        DEFS...,
        label="Arbritrary Units",
        kwargs...
            )

    plot!(plt; xlims=(1e-2, mx), legend=:best, bottom_margin=5mm, kwargs...)
    savefig(plt, fn)
    plt
end

#figure 1.2
figure1_2(
    ;strain="KT2440",
    backbone="pSeva251",
    iptg=2000,
    fn="/home/lewis/thesis/figures/figure1.2.svg",
    legend=:topleft,
    xscale=:log10
);
