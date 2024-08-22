module Hippocampus
using Makie
using Glob
using DrWatson
using DataProcessingHierarchyTools
const DPHT = DataProcessingHierarchyTools
include("utils.jl")
include("eyelink.jl")
include("rawdata.jl")
include("unity.jl")
include("test.jl")
include("sprites.jl")
include("replay.jl")
include("mountainsort.jl")
end