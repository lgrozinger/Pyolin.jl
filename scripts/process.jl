using Pyolin
using DataFrames
using CSV
using Plots
using Plots.Measures


const DATADIR = "/home/campus.ncl.ac.uk/b8051106/sauce/julia/Pyolin.jl/data/"
const FIGDIR = "/home/campus.ncl.ac.uk/b8051106/sauce/julia/Pyolin.jl/figures/"

PerExperiment(outputfn) = Experiments(Pyolin.index) |> DataFrame |> CSV.write(DATADIR * outputfn)

function PerGate(inputfn, outputfn)
    frame = CSV.read(DATADIR * inputfn, DataFrame)
    types = unique(frame[:, [:strain, :backbone, :plasmid]])
    filter!(r -> !(r.plasmid in [Pyolin.INPUT, Pyolin.AUTOFLUOR, Pyolin.STANDARD]), types)
    results = ResponseFunction[]
    for row in eachrow(types)
        inputs  = Experiments(row.strain, row.backbone, Pyolin.INPUT, DATADIR * inputfn)
        outputs = Experiments(row.strain, row.backbone, row.plasmid, DATADIR * inputfn)
        push!(results, ResponseFunction(inputs, outputs))
    end
    results |> DataFrame |> CSV.write(DATADIR * outputfn)
end

function PerGateRpu(inputfn, outputfn)
    frame = CSV.read(DATADIR * inputfn, DataFrame)
    types = unique(frame[:, [:strain, :backbone, :plasmid]])
    filter!(r -> !(r.plasmid in [Pyolin.INPUT, Pyolin.AUTOFLUOR, Pyolin.STANDARD]), types)

    results = ResponseFunction[]
    for row in eachrow(types)
        autos      = Experiments(row.strain, row.backbone, Pyolin.AUTOFLUOR, DATADIR * inputfn)
        standards  = Experiments(row.strain, row.backbone, Pyolin.STANDARD, DATADIR * inputfn)
        inputs     = Experiments(row.strain, row.backbone, Pyolin.INPUT, DATADIR * inputfn)
        outputs    = Experiments(row.strain, row.backbone, row.plasmid, DATADIR * inputfn)
        push!(results, ResponseFunction(
            rpuconvert(inputs, autos, standards),
            rpuconvert(outputs, autos, standards)
        ))
    end
    results |> DataFrame |> CSV.write(DATADIR * outputfn)
end

function Compatibilities(inputfn, outputfn)
    df = CSV.read(DATADIR * inputfn, DataFrame)
    rfs = ResponseFunctions(DATADIR * inputfn)
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

PerExperiment("experimentstats.csv")
PerGate("experimentstats.csv", "gatestats.csv")
PerGateRpu("experimentstats.csv", "gatestatsrpu.csv")
Compatibilities("gatestats.csv", "compatibilities.csv")
Compatibilities("gatestatsrpu.csv", "compatibilitiesrpu.csv")

function compatibilitymap(fn, strain, backbone)
    df = CSV.read(DATADIR * fn, DataFrame)
    df = Pyolin.context_groupings(df; strain=strain, backbone=backbone)

    x, y, z = Pyolin.compatibilitymatrix(df)
    x = (s -> replace(s, "_" => " ")).(x)
    y = (s -> replace(s, "_" => " ")).(y)

    title = "Host: $(strain === missing ? string(:Any) : strain) Backbone: $(backbone === missing ? string(:Any) : backbone)"
    opts = Dict(
        :xrotation => 90,
        :guidefontsize => 12,
        :legendfontsize => 12,
        :tickfontsize => 10,
        :xlabel => "Gate A",
        :ylabel => "Gate B",
        :title => title,
        :bottommargin => 5mm
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

# compatibilitymaps("compatibilities.csv")
# compatibilitymaps("compatibilitiesrpu.csv")
