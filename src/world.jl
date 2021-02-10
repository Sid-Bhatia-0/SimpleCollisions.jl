struct World
    bodies::Vector{RigidBody}
end

get_bodies(world::World) = world.bodies

@pretty_print World

function step!(world::World, dt)
    bodies = get_bodies(world)
    num_bodies = length(bodies)

    # accumulate changes caused by pair-wise interaction of bodies
    for i in 1:num_bodies-1
        for j in i+1:num_bodies
            a = bodies[i]
            b = bodies[j]
            if is_colliding(a, b)
                manifold_ba = Manifold(a, b)
                velocity_change_a, velocity_change_b = resolve_collision(a, b, manifold_ba)
                add_velocity_change!(a, velocity_change_a)
                add_velocity_change!(b, velocity_change_b)
            end
        end
    end

    # update bodies
    for body in bodies
        # apply velocity change
        inv_mass = get_inv_mass(body)
        force = get_force(body)
        add_velocity_change!(body, inv_mass * force * dt)
        apply_velocity_change!(body)

        # apply angular velocity change
        inv_inertia = get_inv_inertia(body)
        torque = get_torque(body)
        add_angular_velocity_change!(body, inv_inertia * torque * dt)
        apply_angular_velocity_change!(body)

        # apply position change
        velocity = get_velocity(body)
        add_position_change!(body, velocity * dt)
        apply_position_change!(body)

        # apply angle change
        angular_velocity = get_angular_velocity(body)
        add_angle_change!(body, angular_velocity * dt)
        apply_angle_change!(body)

        # update axes using current angle
        angle = get_angle(body)
        x_cap = get_x_cap(body)
        new_x_cap = typeof(x_cap)(cos(angle), sin(angle))
        set_x_cap!(body, new_x_cap)
        set_y_cap!(body, rotate_90(new_x_cap))
    end

    return world
end

function simulate!(world::World, num_iter, dt)
    for i in 1:num_iter
        step!(world, dt)
    end

    return world
end
