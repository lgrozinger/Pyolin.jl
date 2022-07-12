using Pyolin
using DataFrames
using ProgressMeter
using CSV

DATADIR = "/home/lewis/sauce/julia/Pyolin/data/"

function ExperimentStatistics(filename)
    df = Pyolin.fileindex[:, [:strain, :backbone, :plasmid, :iptg]]

    es = Vector{Experiment}(undef, length(eachrow(df)))
    @showprogress for (i, r) in enumerate(eachrow(df))
        es[i] = Experiment(r.strain, r.backbone, r.plasmid, r.iptg)
    end
    es |> DataFrame |> CSV.write(DATADIR * filename)
end

function GateStatistics(filename)
    df = Pyolin.fileindex[:, [:strain, :backbone, :plasmid]]
    filter!(r -> r.plasmid ∉ ["1201", "1717", "1818"], df)

    gs = Vector{InputOutput}(undef, length(eachrow(df)))
    @showprogress for (i, r) in enumerate(eachrow(df))
        n = Network(GeneticInverter, r.strain, r.backbone, r.plasmid)
        gs[i] = InputOutput(n)
    end
    gs |> DataFrame |> CSV.write(DATADIR * filename)
end

function GateStatisticsRpu(filename)
    df = Pyolin.fileindex[:, [:strain, :backbone, :plasmid]]
    filter!(r -> r.plasmid ∉ ["1201", "1717", "1818"], df)

    gs = Vector{InputOutput}(undef, length(eachrow(df)))
    @showprogress for (i, r) in enumerate(eachrow(df))
        n = Network(GeneticInverter, r.strain, r.backbone, r.plasmid)
        a = Network(Autofluorescent, r.strain, r.backbone)
        s = Network(Standardisation, r.strain, r.backbone)
        gs[i] = InputOutput(n, a, s)
    end
    gs |> DataFrame |> CSV.write(DATADIR * filename)
end

ExperimentStatistics("experimentstats.csv")
GateStatistics("gatestats.csv")
GateStatisticsRpu("gatestatsrpu.csv")


#GateStats() |> CSV.write("/home/lewis/sauce/julia/Pyolin/data/gatestats.csv")
#statsframerpu |> CSV.write("/home/lewis/sauce/julia/Pyolin/data/statsrpu.csv")

#statsframe = CSV.read("/home/lewis/sauce/julia/Pyolin/data/stats.csv", DataFrame)
#statsframerpu = CSV.read("/home/lewis/sauce/julia/Pyolin/data/statsrpu.csv", DataFrame)

