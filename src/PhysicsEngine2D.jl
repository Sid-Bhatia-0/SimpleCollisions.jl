module PhysicsEngine2D

export PE2D
const PE2D = PhysicsEngine2D

import GeometryBasics
const GB = GeometryBasics
import LinearAlgebra
const LA = LinearAlgebra
import Requires
import MacroTools

include("printing.jl")
include("collision_detection.jl")
include("shapes.jl")
include("rigid_body.jl")
include("manifold.jl")
include("world.jl")

function __init__()
    Requires.@require Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" include("rendering.jl")
end

end
