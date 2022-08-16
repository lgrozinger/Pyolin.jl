using Pyolin
using DataFrames
using CSV
using Plots


DATADIR = "data/"
FIGDIR  = "figures/"

PerExperiment(fn) = Experiments(Pyolin.index) |> DataFrame |> CSV.write(pwd()*"/"*DATADIR*fn)

function PerGateRpu(inputfn, outputfn)
    frame = CSV.read(DATADIR * inputfn, DataFrame)
    types = unique(frame[:, [:strain, :backbone, :plasmid]])
    filter!(r -> r.plasmid ∉ [Pyolin.INPUT, Pyolin.AUTOFLUOR, Pyolin.STANDARD], types)

    results = Hill[]
    for row in eachrow(types)
        inputs  = Experiments(search(row.strain, row.backbone, Pyolin.INPUT))
        autos   = Auto.(inputs)
        stds    = Standard.(inputs)
        outputs = Experiments(search(row.strain, row.backbone, row.plasmid))
        px(mm) = mm * 96 / 25.4
        plt = notplot(
            rpuconvert(inputs, autos, stds),
            rpuconvert(outputs, autos, stds),
            thresholds=true,
            guidefontsize=9,
            legendfontsize=9,
            tickfontsize=9,
            size=(px(135), px(70)),
        )
        savefig(plt, pwd()*"/"*FIGDIR*"response/$(row.strain)-$(row.backbone)-$(row.plasmid)-rpu.svg")
        push!(results, Hill(rpuconvert(inputs, autos, stds), rpuconvert(outputs, autos, stds)))
    end
    results |> DataFrame |> CSV.write(pwd()*"/"*DATADIR*outputfn)
end


PerExperiment("experimentstats.csv")
PerGateRpu("experimentstats.csv", "gatestatsrpu.csv")

