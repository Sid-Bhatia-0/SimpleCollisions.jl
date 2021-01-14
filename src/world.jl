mutable struct World
    bodies::Vector{RigidBody}
end

function update!(world::World, dt)
    for body in world.bodies
        update!(body, dt)
    end
end
