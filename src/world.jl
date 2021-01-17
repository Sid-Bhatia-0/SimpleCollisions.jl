struct World
    bodies::Vector{RigidBody}
end

get_bodies(world::World) = world.bodies

@pretty_print World

function step!(world::World, dt)
    bodies = get_bodies(world)
    num_bodies = length(bodies)

    # updates caused by pair-wise interaction of rigid bodies
    for i in 1:num_bodies-1
        for j in i+1:num_bodies
            a = bodies[i]
            shape_a = get_shape(a)
            b = bodies[j]
            shape_b = get_shape(b)
            if is_colliding(shape_a, shape_b)
                manifold_ba = Manifold(shape_a, shape_b)
                velocity_change_a, velocity_change_b = resolve_collision!(a, b, manifold_ba)
                add_velocity_change!(a, velocity_change_a)
                apply_velocity_change!(a)
                add_velocity_change!(b, velocity_change_b)
                apply_velocity_change!(b)
            end
        end
    end

    # updates caused by self (individual body)
    for body in bodies
        # accumulate velocity change caused by net force
        inv_mass = get_inv_mass(body)
        force = get_force(body)
        add_velocity_change!(body, inv_mass * force * dt)
        apply_velocity_change!(body)

        # accumulate position change caused by net velocity
        velocity = get_velocity(body)
        add_position_change!(body, velocity * dt)
        apply_position_change!(body)

        # update shape caused by change in position
        shape = get_shape(body)
        position = get_position(body)
        set_shape!(body, translated(shape, position))
    end

    return world
end

function simulate!(world::World, num_iter, dt)
    for i in 1:num_iter
        step!(world, dt)
    end

    return world
end
