mutable struct World
    bodies::Vector{RigidBody}
end

get_bodies(world::World) = world.bodies

@pretty_print World

function update!(world::World, dt)
    for body in world.bodies
        update!(body, dt)
    end
    return world
end

function Base.run(world::World, num_iter, dt)
    for i in 1:num_iter
        update!(world, dt)
    end
    return world
end
