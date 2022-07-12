using Pyolin
using Plots
using DataFrames
using CSV
using ProgressMeter
using Combinatorics


function full_matrix(rpu=false)
    compatibilities = DataFrame(
        strain_a = String15[],
        strain_b = String15[],
        backbone_a = String15[],
        backbone_b = String15[],
        plasmid_a = String15[],
        plasmid_b = String15[],
        compatible = Bool[]
    )

    gates = filter(x -> x[:plasmid] ∉ ["1818", "1201", "1717"], Pyolin.ALL)
    n = length(gates) ^ 2
    i = 0
    for a in gates
        gateA = Gate(;a..., rpu=rpu)
        for b in gates
            print("$i / $n\r")
            i = i + 1
            compat = false
            if split(a[:plasmid], "_")[1] != split(b[:plasmid], "_")[1]
                gateB = Gate(;b..., rpu=rpu)
                compat = compatible(gateA, gateB)
            end
            append!(
                compatibilities,
                DataFrame(
                    strain_a = a[:strain],
                    strain_b = b[:strain],
                    backbone_a = a[:backbone],
                    backbone_b = b[:backbone],
                    plasmid_a = a[:plasmid],
                    plasmid_b = b[:plasmid],
                    compatible = [compat]
            ))
        end
    end
    compatibilities
end

C = full_matrix()
C |> CSV.write("/home/lewis/sauce/julia/Pyolin/data/compatibilities.csv")

D = full_matrix(rpu=true)
D |> CSV.write("/home/lewis/sauce/julia/Pyolin/data/compatibilities_rpu.csv")
