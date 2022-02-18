# dynamic body vs. dynamic body
function get_normal_impulse(inv_mass_a, inv_mass_b, initial_velocity_ao::Vector2D, initial_velocity_bo::Vector2D, e, normal_o::Vector2D)
    initial_velocity_ba = initial_velocity_bo .- initial_velocity_ao
    initial_velocity_ba_normal_o = LA.dot(initial_velocity_ba, normal_o)
    final_velocity_ba_normal_o = -e * initial_velocity_ba_normal_o
    final_velocity_ba_normal_o = max(zero(final_velocity_ba_normal_o), final_velocity_ba_normal_o)
    j_ao_normal_o = (initial_velocity_ba_normal_o - final_velocity_ba_normal_o) / (inv_mass_a + inv_mass_b)
    j_bo_normal_o = -j_ao_normal_o
    return j_ao_normal_o, j_bo_normal_o
end

# dynamic body vs. kinematic body
get_normal_impulse(inv_mass_b, initial_velocity_ao::Vector2D, initial_velocity_bo::Vector2D, e, normal_o::Vector2D) = get_normal_impulse(zero(inv_mass_b), inv_mass_b, initial_velocity_ao, initial_velocity_bo, e, normal_o)

# dynamic body vs. static body
get_normal_impulse(inv_mass_b, initial_velocity_bo::Vector2D, e, normal_o::Vector2D) = get_normal_impulse(zero(inv_mass_b), inv_mass_b, zero(Vector2D), initial_velocity_bo, e, normal_o)

function get_tangential_impulse(inv_mass_a, inv_mass_b, initial_velocity_ao::Vector2D, initial_velocity_ba::Vector2D, tangent_o::Vector2D, mu_s, mu_k, j_ao_normal_o)
    initial_velocity_ba = initial_velocity_bo .- initial_velocity_ao
    initial_velocity_ba_tangent_o = LA.dot(initial_velocity_ba, tangent_o)

    max_j_ao_tangent_o = initial_velocity_ba_tangent_o / (inv_mass_a + inv_mass_b)

    if abs(max_j_ao_tangent_o) <= abs(mu_s * j_ao_normal_o)
        j_ao_tangent_o = max_j_ao_tangent_o
        j_bo_tangent_o = -j_ao_tangent_o
        return j_ao_tangent_o, j_bo_tangent_o
    else
        j_ao_tangent_o = sign(max_j_ao_tangent_o) * abs(mu_k * j_ao_normal_o)
        j_bo_tangent_o = -j_ao_tangent_o
        return j_ao_tangent_o, j_bo_tangent_o
    end
end
