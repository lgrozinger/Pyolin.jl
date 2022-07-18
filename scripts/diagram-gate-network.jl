using Pyolin
using Plots
using StatsPlots
using DataFrames

const FIGDIR = "/home/lewis/sauce/julia/Pyolin/figures/"

function gatecytometry()
    justgates = filter(r -> !(r.plasmid ∈ ["1201", "1717", "1818"]), Pyolin.index)
    gates = groupby(justgates, [:strain, :backbone, :plasmid])

    px(mm) = mm * 96 / 25.4

    for gate in gates
        filtered = filter(r -> r.iptg ∈ [0, 5, 10, 50, 100, 500, 1000], gate)
        experiments = Experiments(filtered)

        plt = plot(
            xlabel="YFP",
            ylabel="Density",
            legendfontsize=10,
            guidefontsize=10,
            dpi=900,
            size=(px(100), px(60)),
            grid=false
        )
        xlim = maximum(quantile.(experiments, 0.95))
        plot!(plt, xlims=(0, xlim), palette=palette(:Accent_8))

        for (i, e) in enumerate(experiments)
            density!(plt, e, label="IPTG=$(e.iptg)", linewidth=3, color=i)
            vline!(plt, [median(e)], color=i, label="", linestyle=:dash)
        end

        strain = experiments[1].strain
        backbone = experiments[1].backbone
        plasmid = experiments[1].plasmid
        savefig(plt, FIGDIR * "experiments/$(strain)-$(backbone)-$(plasmid)-densities.svg")
    end
end

function inputcytometry()
    justgates = filter(r -> r.plasmid == "1818", Pyolin.index)
    gates = groupby(justgates, [:strain, :backbone])

    px(mm) = mm * 96 / 25.4

    for gate in gates
        filtered = filter(r -> r.iptg ∈ [0, 5, 10, 50, 100, 500, 1000], gate)
        experiments = Experiments(filtered)

        plt = plot(
            xlabel="YFP",
            ylabel="Density",
            legendfontsize=10,
            guidefontsize=10,
            dpi=900,
            size=(px(100), px(60)),
            grid=false
        )
        xlim = maximum(quantile.(experiments, 0.95))
        plot!(plt, xlims=(0, xlim), palette=palette(:Accent_8))

        for (i, e) in enumerate(experiments)
            density!(plt, e, label="IPTG=$(e.iptg)", linewidth=3, color=i)
            vline!(plt, [median(e)], color=i, label="", linestyle=:dash)
        end

        strain = experiments[1].strain
        backbone = experiments[1].backbone
        plasmid = experiments[1].plasmid
        savefig(plt, FIGDIR * "experiments/$(strain)-$(backbone)-$(plasmid)-densities.svg")
    end
end
