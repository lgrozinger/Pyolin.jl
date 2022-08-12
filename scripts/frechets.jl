using Pyolin
using DataFrames
using CSV
using Plots
using Plots.Measures

DATADIR = "/home/lewis/sauce/julia/Pyolin/data/"
FIGDIR = "/home/lewis/sauce/julia/Pyolin/figures/"


function polygonal_curve(fn, row)
    inputs = Experiments(row.strain, row.backbone, Pyolin.INPUT, fn)
    outputs = Experiments(row.strain, row.backbone, row.plasmid, fn)
    # autos = Autos(inputs)
    standards = Experiments(row.strain, row.backbone, Pyolin.STANDARD, fn)

    x = median.(rpuconvert(inputs, standards))
    y = median.(rpuconvert(outputs, standards))
    hcat(x, y)'
end

function frechet_distance(fn, rowA, rowB)
    A = polygonal_curve(fn, rowA)
    B = polygonal_curve(fn, rowB)
    Pyolin.frechet(A, B)
end

function frechet_matrix(fn, frame)
    N, _ = size(frame)
    M = Matrix{Float64}(undef, N, N)

    for (j, a) in enumerate(eachrow(frame))
        for (i, b) in enumerate(eachrow(frame))
            M[i, j] = frechet_distance(fn, a, b)
        end
    end
    M
end

function frechet_frame(frame, M)
    strains = frame.strain
    backbones = frame.backbone

    n, _ = size(M)

    DataFrame(
        strainA = repeat(strains, n),
        strainB = reduce(vcat, repeat([strains[i]], n) for i in 1:n),
        backboneA = repeat(backbones, n),
        backboneB = reduce(vcat, repeat([backbones[i]], n) for i in 1:n),
        distance = M[:]
    )
end

function frechet_map(fn, frames...)
    frames = [sort(frame, [:strain, :backbone]) for frame in frames]
    zs = [frechet_matrix(fn, frame) for frame in frames]
    z = reduce(.+, zs)
    z = z ./ maximum(z)
    n,m = size(z)

    x = frames[1].strain .* "\n" .* frames[1].backbone

    px(mm) = mm * 96 / 25.4

    opts = Dict(
        :xrotation => 90,
        :guidefontsize => 9,
        :tickfontsize => 9,
        :aspectratio => 1,
        :cbar => false,
        :xlims => (0, size(z)[1]),
        :ylims => (0, size(z)[2]),
        :size => (px(100), px(100)),
        :showaxes => false,
    )
    plt = plot()
    heatmap!(
        plt,
        collect(0:length(x)),
        collect(0:length(x)),
        z;
        xticks=(collect(0:length(x)) .+ 0.5, x),
        yticks=(collect(0:length(x)) .+ 0.5, x),
        opts...)

    anots = [text(round(m; digits=2), "Computer Modern"; pointsize=9, color=(m < maximum(z) / 2 ? :white : :black)) for m in z[:]]
    scatter!(
        plt,
        repeat(1:n, inner=n) .- 0.5, repeat(1:m, outer=m) .- 0.5;
        markerstrokecolor = RGBA(0, 0, 0, 0.0),
        seriesalpha = 0.0,
        label="",
        series_annotations = anots,
        opts...
    )
    plt
end

gatefile = DATADIR * "gatestats.csv"
exprfile = DATADIR * "experimentstats.csv"
gates = groupby(CSV.read(gatefile, DataFrame), [:plasmid])
gates = [gate for gate in gates if size(gate) == (7, 7)]

for gate in gates
    plt = frechet_map(exprfile, gate)
    savefig(plt, FIGDIR * "$(gate.plasmid[1])-frechets.svg")
end

gates = groupby(CSV.read(gatefile, DataFrame), [:plasmid])
gates = [gate for gate in gates if size(gate) == (7, 7)]
plt = frechet_map(exprfile, gates...)
savefig(plt, FIGDIR * "summed-frechets.svg")
