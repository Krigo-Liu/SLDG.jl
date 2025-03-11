module SLDG

# greet() = print("Hello World!")
include("SLDG1d.jl")
using .SLDG1d

export SimulationParameters, SimulationState 
export Vertex, EulerianElement, Segment, UpstreamElement
export sldg1d

include("../test/test1d.jl")


end # module SLDG
