module PhysicsPrimitives2D

import StaticArrays
const SA = StaticArrays
import LinearAlgebra
const LA = LinearAlgebra

include("utils.jl")
include("position_and_orientation.jl")
include("shapes.jl")
include("collision_detection.jl")
include("collision_manifold.jl")
include("physics.jl")

end
