module Pyolin

using Base.Iterators
using InlineStrings
using CSV
using Tables
using DataFrames
using FCSFiles
using FileIO
using RecipesBase
using Distributions
using Turing
using MCMCChains
using RecipesBase
using StatsPlots
using StatsBase
using LsqFit
using Logging
using Symbolics
using CategoricalArrays

import StatsBase.mean
import StatsBase.median
import StatsBase.var
import StatsBase.minimum
import StatsBase.maximum

import Distributions.fit

import Tables.istable
import Tables.rowaccess
import Tables.rows
import Tables.schema
import Tables.getcolumn
import Tables.columnnames

# configurations
const DATADIR   = "/home/lewis/data/pyolin_dataset/"
const CHANNEL   = "B1-H"
const AUTOFLUOR = "1201"
const STANDARD  = "1717"
const INPUT     = "1818"

include("utils.jl")
include("processing.jl")
export search, events

include("gates.jl")

include("experiment.jl")
export Experiment, Experiments
export median, mean, var, fit

include("responsefunction.jl")
export ResponseFunction, ResponseFunctions
export compatible, rpuconvert

include("model.jl")
export model1717, model1818, cond1717, cond1818
export StandardModel, InputModel

include("plots.jl")

end # module
