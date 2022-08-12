using Pyolin
using Plots
using StatsPlots
using DataFrames

FIGDIR = "/home/lewis/sauce/julia/Pyolin/figures/experiments/"

px(mm) = mm * 96 / 25.4

function cytometry_densities_plot(dataframe)
    experiments = Experiments(dataframe)
    plt = plot(
        xlabel="YFP",
        ylabel="Density",
        legendfontsize=9,
        guidefontsize=9,
        tickfontsize=9,
        dpi=900,
        size=(px(86.5), px(60.0)),
        grid=false,
    )

    for (i, e) in enumerate(experiments)
        density!(plt, e, label="$(e.iptg)" * L"$\mu$M", linewidth=4, color=i)
        vline!(plt, [median(e)], label="", linestyle=:dash, linewidth=2, color=i)
    end

    xlim = maximum(quantile.(experiments, 0.95))
    plot!(plt, xlims=(0, xlim))
    fn = "$(experiments[1].strain)-$(experiments[1].backbone)-$(experiments[1].plasmid)-densities.svg"
    savefig(plt, FIGDIR * fn)
    plt
end

function gatecytometry()
    justgates = filter(r -> !(r.plasmid ∈ ["1201", "1717", "1818"]), Pyolin.index)
    gates = groupby(justgates, [:strain, :backbone, :plasmid])

    for gate in gates
        filter(r -> r.iptg ∈ [0, 5, 10, 50, 100, 500, 1000], gate)
        cytometry_densities_plot(filter(r -> r.iptg ∈ [0, 5, 10, 50, 100, 500, 1000], gate))
    end
end

function inputcytometry()
    just = filter(r -> r.plasmid == "1818", Pyolin.index)
    gates = groupby(just, [:strain, :backbone])

    px(mm) = mm * 96 / 25.4

    for gate in gates
        filter(r -> r.iptg ∈ [0, 5, 10, 50, 100, 500, 1000], gate)
        cytometry_densities_plot(filter(r -> r.iptg ∈ [0, 5, 10, 50, 100, 500, 1000], gate))
    end
end

function autocytometry()
    just = filter(r -> r.plasmid == "1201", Pyolin.index)
    gates = groupby(just, [:strain, :backbone])

    px(mm) = mm * 96 / 25.4

    for gate in gates
        filtered = filter(r -> r.iptg ∈ [0, 5, 10, 50, 100, 500, 1000], gate)
        cytometry_densities_plot(filter(r -> r.iptg ∈ [0, 5, 10, 50, 100, 500, 1000], gate))
    end
end

function stdcytometry()
    just = filter(r -> r.plasmid == "1717", Pyolin.index)
    gates = groupby(just, [:strain, :backbone])

    px(mm) = mm * 96 / 25.4

    for gate in gates
        filtered = filter(r -> r.iptg ∈ [0, 5, 10, 50, 100, 500, 1000], gate)
        cytometry_densities_plot(filter(r -> r.iptg ∈ [0, 5, 10, 50, 100, 500, 1000], gate))
    end
end
