using Pyolin
using Plots


function input_plasmid_plot(;kwargs...)
    plt = plot(;xlabel="Signal (A.U.)", ylabel="Probability", kwargs...)
    for i in Pyolin.IPTGS
        E = Experiment(;strain="KT2440", backbone="pSeva221", plasmid="1818", iptg=i, rpu=false)
        plot!(plt, E.distribution, label="IPTG=$i")
    end
    plt
end


