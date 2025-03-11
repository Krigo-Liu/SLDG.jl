module SLDG

# greet() = print("Hello World!")
include("SLDG1d.jl")
include("test1d.jl")

using .SLDG1d

export SimulationParameters, SimulationState 
export Vertex, EulerianElement, Segment, UpstreamElement
export sldg1d

# 动态加载 test1d.jl
const test1d = let file = joinpath(dirname(@__FILE__), "..", "test", "test1d.jl")
    eval(Meta.parse(read(file, String)))
end

end # module SLDG
