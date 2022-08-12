using Pyolin
using DataFrames
using CSV
using Plots
using StatsPlots
using LinearAlgebra
using Distributions
using LsqFit
using Symbolics
using LaTeXStrings

DATADIR = "/home/lewis/sauce/julia/Pyolin/data/"
FIGDIR = "/home/lewis/sauce/julia/Pyolin/figures/compare-hill-functions/"

px(mm) = mm * 96 / 25.4

frame = CSV.read(DATADIR * "experimentstats.csv", DataFrame)
a = Hill("KT2440", "pSeva221", "Srpr_s1")
b = Hill("KT2440", "pSeva221", "Qacr_q1")

plta = plot(guidefontsize=9, tickfontsize=9, legendfontsize=9, grid=false, palette=:Accent_8, ylims=(0, 1), xlabel="Input", ylabel="Output", size=(px(70), px(60)))
xs = range(0, 1, 100)
plot!(plta, xs, a.(xs), label="Hill model", linewidth=3)
plot!(plta, [Pyolin.outputhigh(a)], seriestype=:hline, label="High Output", linestyle=:dash, linewidth=3)
plot!(plta, [Pyolin.outputlow(a)], seriestype=:hline, label="Low Output", linestyle=:dash, linewidth=3)

pltb = plot(guidefontsize=9, tickfontsize=9, legendfontsize=9, grid=false, palette=:Accent_8, ylims=(0, 1), xlabel="Output", ylabel="Input", size=(px(70), px(60)))
xs = range(0, 1, 100)
plot!(pltb, b.(xs), xs, label="Hill model", linewidth=3)
plot!(pltb, [Pyolin.inputhigh(b)], seriestype=:hline, label="High Input", linestyle=:dash, linewidth=3)
plot!(pltb, [Pyolin.inputlow(b)], seriestype=:hline, label="Low Input", linestyle=:dash, linewidth=3)

savefig(plta, FIGDIR * "inputgate.svg")
savefig(pltb, FIGDIR * "outputgate.svg")







