module PhysicsEngine2D

export PE2D
const PE2D = PhysicsEngine2D

import GeometryBasics
const GB = GeometryBasics
import LinearAlgebra
const LA = LinearAlgebra

include("collision_detection.jl")
include("shapes.jl")

end
