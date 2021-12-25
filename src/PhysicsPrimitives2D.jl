module PhysicsPrimitives2D

import LinearAlgebra as LA
import StaticArrays as SA

include("utils.jl")
include("position_and_orientation.jl")
include("shapes.jl")
include("collision_detection.jl")
include("collision_manifold.jl")
include("physics.jl")

end
