module Pyolin

using InlineStrings
using CSV
using Tables
using DataFrames
using FCSFiles
using FileIO
using RecipesBase
using Distributions
using LaTeXStrings
using Plots
using HypothesisTests

using Turing
using Zygote
Turing.setadbackend(:forwarddiff)

using MCMCChains
using StatsBase
using LsqFit
using Symbolics

import StatsBase.mean
import StatsBase.median
import StatsBase.var
import StatsBase.minimum
import StatsBase.maximum
import StatsBase.sample

import Distributions.fit

import Tables.istable
import Tables.rowaccess
import Tables.rows
import Tables.schema
import Tables.getcolumn
import Tables.columnnames

# configurations
const DATADIR   = "/home/lewis/data/pyolin_dataset/"
const CHANNEL   = "B1-A"
const AUTOFLUOR = "1201"
const STANDARD  = "1717"
const INPUT     = "1818"

include("utils.jl")
include("processing.jl")
export search, events

include("gates.jl")

include("experiment.jl")
export Experiment, Experiments
export median, mean, var, quantile, fit, sample

include("responsefunction.jl")
export ResponseFunction, ResponseFunctions
export compatible, rpuconvert

include("model.jl")
export Standardisation, InputSensor


include("plots.jl")

end # module
