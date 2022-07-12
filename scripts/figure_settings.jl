using Plots
using ColorSchemes

const CSVDIR = "/home/lewis/sauce/julia/Pyolin/data/"
const GPDIR = "/home/lewis/sauce/julia/Pyolin/gnuplots/"

const PALETTE = palette(:Accent_8)

const DEFS = Dict(
    :guidefontsize => 12,
    :legendfontsize => 10,
    :tickfontsize => 10,
    :palette => PALETTE
)

