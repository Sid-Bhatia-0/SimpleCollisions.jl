function resolve_collision!(a::RigidBody, b::RigidBody, manifold_ba::Manifold)
    velocity_a = get_velocity(a)
    velocity_b = get_velocity(b)
    velocity_ba = velocity_b .- velocity_a

    inv_mass_a = get_inv_mass(a)
    inv_mass_b = get_inv_mass(b)

    T = eltype(velocity_a)

    velocity_change_a = zero(velocity_a)
    velocity_change_b = zero(velocity_b)

    normal_ba = get_normal(manifold_ba)
    velocity_ba_normal_ba = LA.dot(velocity_ba, normal_ba)

    if velocity_ba_normal_ba < zero(T)
        # linear impulse
        e = min(get_restitution(a), get_restitution(b))
        # j is the magnitude of impulse. j >= 0
        j = - ((1 + e) * velocity_ba_normal_ba) / (inv_mass_a + inv_mass_b)
        velocity_change_a = velocity_change_a .- j * normal_ba * inv_mass_a
        velocity_change_b = velocity_change_b .+ j * normal_ba * inv_mass_b

        # angular impulse
        tangent_ba = typeof(normal_ba)(-normal_ba[2], normal_ba[1])
        velocity_ba_tangent_ba = LA.dot(velocity_ba, tangent_ba)
        # tangent_ba should be such that velocity_ba_tangent_ba is >= 0
        if velocity_ba_tangent_ba < zero(T)
            tangent_ba = -tangent_ba
            velocity_ba_tangent_ba = -velocity_ba_tangent_ba
        end

        jt = velocity_ba_tangent_ba / (inv_mass_a + inv_mass_b) # jt will be >= 0
        static_friction_coeff = (get_static_friction_coeff(a) + get_static_friction_coeff(b)) / 2
        if abs(jt) < static_friction_coeff * j
            velocity_change_a = velocity_change_a .+ jt * tangent_ba * inv_mass_a
            velocity_change_b = velocity_change_b .- jt * tangent_ba * inv_mass_b
        else
            dynamic_friction_coeff = (get_dynamic_friction_coeff(a) + get_dynamic_friction_coeff(b)) / 2
            velocity_change_a = velocity_change_a .+ j * tangent_ba * inv_mass_a * dynamic_friction_coeff
            velocity_change_b = velocity_change_b .- j * tangent_ba * inv_mass_b * dynamic_friction_coeff
        end
    end

    return (velocity_change_a, velocity_change_b)
end
