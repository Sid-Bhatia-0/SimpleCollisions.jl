module PhysicsPrimitives2D

export PP2D
const PP2D = PhysicsPrimitives2D

import StaticArrays
const SA = StaticArrays
import LinearAlgebra
const LA = LinearAlgebra

include("utils.jl")
include("axes.jl")
include("shapes.jl")
include("collision_detection.jl")
include("collision_manifold.jl")
include("abstract_body.jl")
include("physics.jl")

end
