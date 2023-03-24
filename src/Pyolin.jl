module Pyolin

using CSV
using Tables
using DataFrames
using FCSFiles
using FileIO
using RecipesBase
using Distributions
using LaTeXStrings
using Plots
using StatsPlots
using HypothesisTests
using LinearAlgebra
using Turing
using DynamicHMC
using Zygote
using ReverseDiff
using Memoization
Turing.setadbackend(:forwarddiff)
using MCMCChains
using MCMCChainsStorage
using HDF5
using StatsBase
using Statistics
using Distances
using LsqFit
using Symbolics
using Flux

import StatsBase.mean
import StatsBase.median
import StatsBase.var

import StatsBase.minimum
import StatsBase.maximum
import StatsBase.sample

import HypothesisTests.ApproximateTwoSampleKSTest

import Distributions.fit

import Tables.istable
import Tables.rowaccess
import Tables.rows
import Tables.schema
import Tables.getcolumn
import Tables.columnnames

# configurations
const DATADIR   = "/home/lewis/data/pyolin_dataset/" # where the dataset can be found
const CHANNEL   = "B1-H" # the interesting channel
const AUTOFLUOR = "1201" # the codename for the autoflourescent measurement construct
const STANDARD  = "1717" # the codename for the standardisation construct
const INPUT     = "1818" # the codename for the input sensor measurement construct

include("utils.jl")
include("processing.jl")
include("gates.jl")

include("experiment.jl")
export events
export Experiment, RawExperiment, RPUExperiment
export Auto, Autos, Standard, Standards, Input, Inputs

include("responsefunction.jl")
export Hill, Hills
export compatible, orthogonal
export inputs, outputs, autos, standards
export inputhigh, inputlow, outputhigh, outputlow

include("model.jl")
export save, load
export experiment, strain, backbone, iptg, plasmid, chain, model
export simulate, priorsimulate, priorpredict, posteriorpredict, conditions

export Constitutive

#include("input.jl")
#export Input

#include("nots.jl")
#export NotGate

include("plots.jl")

export median, mean, var, cov, quantile, fit, sample, varvar, ApproximateTwoSampleKSTest

# include("probabilistic.jl")

end # module
