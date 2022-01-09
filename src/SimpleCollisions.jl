module SimpleCollisions

import LinearAlgebra as LA
import StaticArrays as SA

const Vector2D{T} = SA.SVector{2, T}

include("utils.jl")
include("position_and_orientation.jl")
include("standard_shapes.jl")
include("collision_detection.jl")
include("collision_manifold.jl")
include("physics.jl")

end
