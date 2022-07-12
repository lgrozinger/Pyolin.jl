using Pyolin
using Distributions
using StatsPlots
using Plots.PlotMeasures


include("figure_settings.jl")

function figure1_3_A(;fn, kwargs...)
    plt = plot(
        0:0.01:0.5,
        ones(length(0:0.01:0.5));
        DEFS...,
        label="Ideal",
        linewidth=16
    )
    plot!(
        plt,
        0.5:0.01:1,
        zeros(length(0.5:0.01:1));
        DEFS...,
        label="",
        color=1,
        linewidth=16
    )
    plot!(
        plt,
        [0.5, 0.5],
        [1, 0];
        DEFS...,
        label="",
        color=1,
        style=:dash,
        linewidth=6
    )
    plot!(
        plt,
        0:0.01:1,
        (x -> 0.95 / ((2 * x)^4 + 1)).(0:0.01:1);
        linewidth=6,
        label="Practical",
        color=3
    )
    plot!(
        plt;
        DEFS...,
        xlims=(0, 1),
        ylims=(0, 1),
        bottom_margin=5mm,
        xlabel="A expression",
        ylabel="Z expression",
        size=(400, 250),
        xticks=[0, 1/2, 1],
        yticks=[0, 1/2, 1],
        kwargs...
            )
    savefig(plt, fn)
    plt
end

function figure1_3_B(;fn, kwargs...)
    g = Gate(strain="KT2440", backbone="pSeva221", plasmid="Srpr_s1", rpu=true)
    plt = plotgate(g; title="", xlabel="Input (RPU)", ylabel="Output (RPU)", DEFS...)
    savefig(plt, fn)
    plt
end

#figure 1.3
figure1_3_A(
    fn="/home/lewis/thesis/figures/figure1.3.A.svg",
);

figure1_3_B(
    fn="/home/lewis/thesis/figures/figure1.3.B.svg",
);
