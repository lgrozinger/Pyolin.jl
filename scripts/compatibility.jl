using Pyolin
using DataFrames
using CSV
using Plots
using Plots.Measures


DATADIR = "/home/lewis/sauce/julia/Pyolin/data/"
FIGDIR = "/home/lewis/sauce/julia/Pyolin/figures/"


function Compatibilities(inputfn, outputfn)
    df = CSV.read(DATADIR * inputfn, DataFrame)
    rfs = Hills(df)
    sA, sB, bA, bB, pA, pB, c = [], [], [], [], [], [], []
    for a in rfs
        for b in rfs
            push!(sA, a.strain)
            push!(sB, b.strain)
            push!(bA, a.backbone)
            push!(bB, b.backbone)
            push!(pA, a.plasmid)
            push!(pB, b.plasmid)
            push!(c, compatible(a, b))
        end
    end
    DataFrame(sA=sA, sB=sB, bA=bA, bB=bB, pA=pA, pB=pB, c=c) |> CSV.write(DATADIR * outputfn)
end

function compatibilitymap(fn, strain, backbone)
    df = CSV.read(DATADIR * fn, DataFrame)
    df = Pyolin.context_groupings(df; strain=strain, backbone=backbone)

    x, y, z = Pyolin.compatibilitymatrix(df)
    x = (s -> replace(s, "_" => " ")).(x)
    y = (s -> replace(s, "_" => " ")).(y)

    title = "Host: $(strain === missing ? string(:Any) : strain) Backbone: $(backbone === missing ? string(:Any) : backbone)"

    println(title)
    println(sum((x -> x == "Compatible" ? 1 : 0).(z)))

    px(mm) = mm * 96 / 25.4
    opts = Dict(
        :xrotation => 90,
        :guidefontsize => 9,
        :legendfontsize => 9,
        :tickfontsize => 9,
        :xlabel => "Gate A",
        :ylabel => "Gate B",
        :size => (px(135), px(120)),
        :aspectratio => 1,
    )
    categorymap(x, y, z; opts...)
end

function compatibilitymaps(fn)
    strains = [missing; unique(Pyolin.index.strain)]
    backbones = [missing; unique(Pyolin.index.backbone)]

    for strain in strains
        for backbone in backbones
            plt = compatibilitymap(fn, strain, backbone)
            savefig(plt, FIGDIR * "compatible/$(strain)-$(backbone)-$(fn).svg")
        end
    end
end

Compatibilities("gatestatsrpu.csv", "compatibilitiesrpu.csv")
compatibilitymaps("compatibilitiesrpu.csv")

