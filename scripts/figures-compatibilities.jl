using Pyolin
using DataFrames
using CSV
using Plots

DATADIR = "/home/lewis/sauce/julia/Pyolin/data/"
FIGDIR = "/home/lewis/sauce/julia/Pyolin/figures/compatible/"

function mapping(df; strain=missing, backbone=missing, kwargs...)
    function formatname(name)
        rep, rbs = split(name, "_")
        rep * " " * rbs
    end
    data = Pyolin._context_groupings(df; strain=strain, backbone=backbone)
    groups = groupby(data, [:plasmidA])

    x = formatname.([g[1, :plasmidA] for g in groups])
    y = x
    z = [[r.compatible for r in eachrow(group)] for group in groups]
    z = reduce(hcat, z)
    compatmap(x, y, z; kwargs...)
end

df = CSV.read(DATADIR * "compatibilities.csv", DataFrame)

for ctx in Pyolin.CONTEXTS
    plt = mapping(df; ctx...)
    savefig(plt, FIGDIR * "$(ctx[:strain])-$(ctx[:backbone])-compat.svg")
end

for strain in Pyolin.STRAINS
    plt = mapping(df; strain=strain, backbone=missing)
    savefig(plt, FIGDIR * "$(strain)-compat.svg")
end

plt = mapping(df; strain=missing, backbone=missing)
savefig(plt, FIGDIR * "all-compat.svg")

df = CSV.read(DATADIR * "compatibilitiesrpu.csv", DataFrame)

for ctx in Pyolin.CONTEXTS
    plt = mapping(df; ctx...)
    savefig(plt, FIGDIR * "$(ctx[:strain])-$(ctx[:backbone])-compat-rpu.svg")
end

for strain in Pyolin.STRAINS
    plt = mapping(df; strain=strain, backbone=missing)
    savefig(plt, FIGDIR * "$(strain)-compat-rpu.svg")
end

plt = mapping(df; strain=missing, backbone=missing)
savefig(plt, FIGDIR * "all-compat-rpu.svg")

