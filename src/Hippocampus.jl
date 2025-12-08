module Hippocampus
using Makie
using Glob
using DrWatson
using DataProcessingHierarchyTools
using RippleTools
using JLD2
using ProgressMeter
using Meshes
using CRC32c

const DPHT = DataProcessingHierarchyTools
include("topology.jl")
include("utils.jl")
include("paths.jl")
include("sprites.jl")
include("neuropixels.jl")
include("rawdata.jl")
include("eyelink.jl")
include("spikedata.jl")
include("unity.jl")
include("spatial.jl")
include("non_spatial.jl")
include("replay.jl")
include("models.jl")
include("decoding.jl")
#include("mountainsort.jl") # until InteractiveViz is compatible with Makie 0.24
include("glmfits.jl")

include("poster_analysis.jl")
end
