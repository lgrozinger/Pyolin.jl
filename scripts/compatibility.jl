using Pyolin
using DataFrames
using ProgressMeter
using CSV

DATADIR = "/home/lewis/sauce/julia/Pyolin/data/"

function _hill(row)
    y0, y1, k, n = row.ymin, row.ymax, row.khalf, row.ncoeff
    x -> y0 + (y1 - y0) * k^n / (x^n + k^n)
end

function _orthogonal(A, B)
    A.strain != B.strain || split(A.plasmid, "_") != split(B.plasmid, "_")
end

function Compatibilities(filename)
    df = CSV.read(DATADIR * "gatestats.csv", DataFrame)
    
    function gatepair(rowA, rowB)
        f = _hill(rowB)
        connect = (
            _orthogonal(rowA, rowB)
            && rowA.outputlow * 2 < rowA.outputhigh
            && rowB.outputlow * 2 < rowB.outputhigh
            && f(rowA.outputhigh) < rowB.outputlow
            && f(rowA.outputlow) > rowB.outputhigh
        )
        DataFrame(
            strainA = [rowA.strain)],
            strainB = [rowB.strain)],
            backboneA = [rowA.backbone],
            backboneB = [rowB.backbone],
            plasmidA = [rowA.plasmid],
            plasmidB = [rowB.plasmid],
            compatible = [connect]
        )
    end

    results = DataFrame()

    @showprogress for rowA in eachrow(df)
        for rowB in eachrow(df)
            append!(results, gatepair(rowA, rowB))
        end
    end
    results |> CSV.write(DATADIR * filename)
end

function CompatibilitiesRpu(filename)
    df = unique(Pyolin.fileindex[:, [:strain, :backbone, :plasmid]])
    filter!(r -> r.plasmid ∉ ["1201", "1717", "1818"], df)

    function gatepair(a, b)
        DataFrame(
            strainA = [strain(a)],
            strainB = [strain(b)],
            backboneA = [backbone(a)],
            backboneB = [backbone(b)],
            plasmidA = [plasmid(a)],
            plasmidB = [plasmid(b)],
            compatible = [compatible(a, b)]
        )
    end

    results = DataFrame()

    @showprogress for (s, b, p) in eachrow(df)
        n = Network(GeneticInverter, s, b, p)
        auto = Network(Autofluorescent, s, b)
        std = Network(Standardisation, s, b)
        a = InputOutput(n, auto, std)
        for (t, c, q) in eachrow(df)
            m = Network(GeneticInverter, t, c, q)
            b = InputOutput(m, auto, std)
            append!(results, gatepair(a, b))
        end
    end
    results |> CSV.write(DATADIR * filename)
end

Compatibilities("compatibilities.csv")
#CompatibilitiesRpu("compatibilitiesrpu.csv")
