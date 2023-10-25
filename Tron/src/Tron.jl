module Tron # Traffic Road On Networks

### Packages
using Reexport
@reexport using LaTeXStrings
@reexport using Random
@reexport using MultiFloats

@reexport using LinearAlgebra
@reexport using ForwardDiff
@reexport using PreallocationTools

@reexport using Graphs: vertices, outneighbors, inneighbors, dst, src, add_edge!, edges, SimpleDiGraph

@reexport using Plots

@reexport using DelimitedFiles
@reexport using Images: @colorant_str

@reexport using Test
@reexport using TickTock
@reexport using Statistics

### Sources

# include("Displays.jl")
include("Tests.jl")
include("Junction.jl")
include("Network.jl")
# include("OptiNetwork.jl")

include("InitialData.jl")
include("Explicit.jl")

include("Model.jl")
include("Diagnostics.jl")
include("Optim.jl")

include("Utils.jl")
include("MyPlots.jl")



end #module 