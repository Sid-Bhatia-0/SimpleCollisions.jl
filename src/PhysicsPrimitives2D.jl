module PhysicsPrimitives2D

export PP2D
const PP2D = PhysicsPrimitives2D

import GeometryBasics
const GB = GeometryBasics
import LinearAlgebra
const LA = LinearAlgebra

include("utils.jl")
include("axes.jl")
include("shapes.jl")
include("collision_detection.jl")
include("manifold.jl")

end
