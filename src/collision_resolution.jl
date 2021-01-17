function resolve_collision!(a::RigidBody, b::RigidBody, manifold_ba::Manifold)
    normal_ba = get_normal(manifold_ba)
    velocity_a = get_velocity(a)
    velocity_b = get_velocity(b)
    velocity_ba = velocity_b .- velocity_a

    velocity_ba_normal_ba = LA.dot(velocity_ba, normal_ba)

    if velocity_ba_normal_ba > zero(eltype(velocity_a))
        return (zero(velocity_a), zero(velocity_b))
    else
        println("Need to resolve collision!")
        e = min(get_restitution(a), get_restitution(b))

        inv_mass_a = get_inv_mass(a)
        inv_mass_b = get_inv_mass(b)

        j = - ((1 + e) * velocity_ba_normal_ba) / (inv_mass_a + inv_mass_b)
        velocity_change_a = - j * normal_ba * inv_mass_a
        velocity_change_b = j * normal_ba * inv_mass_b

        return (velocity_change_a, velocity_change_b)
    end
end
