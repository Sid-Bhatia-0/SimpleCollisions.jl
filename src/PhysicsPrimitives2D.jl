module PhysicsEngine2D

export PE2D
const PE2D = PhysicsEngine2D

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
