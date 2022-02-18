struct Manifold{T}
    penetration::T
    normal::Vector2D
    contact::Vector2D
end

get_penetration(manifold::Manifold) = manifold.penetration
get_normal(manifold::Manifold) = manifold.normal
get_tangent(manifold::Manifold) = rotate_minus_90(get_normal(manifold))
get_contact(manifold::Manifold) = manifold.contact

#####
# StandardCircle vs. StandardCircle
#####

function Manifold(a::StandardCircle, b::StandardCircle, pos_ba)
    r_a = get_radius(a)
    r_b = get_radius(b)
    T = typeof(r_a)
    normal = LA.normalize(pos_ba)
    if all(isfinite.(normal))
        penetration = r_a + r_b - LA.norm(pos_ba)
        contact = (r_a - penetration / 2) * normal
        return Manifold(penetration, normal, contact)
    else
        penetration = r_a + r_b
        normal = Vector2D(one(T), zero(T))
        contact = (r_a - penetration / 2) * normal
        return Manifold(penetration, normal, contact)
    end
end

Manifold(a::StandardCircle, b::StandardCircle, pos_ba, dir_ba) = Manifold(a, b, pos_ba)

#####
# HyperRectangle vs. HyperSphere
#####

function get_closest_edge_from_inside(rect::StandardRect, pos)
    half_width = get_half_width(rect)
    half_height = get_half_height(rect)

    x = pos[1]
    y = pos[2]

    penetration, edge_id = findmin((y + half_height, half_width - x, half_height - y, x + half_width))

    return (penetration, edge_id)
end

function Manifold(a::StandardRect, b::StandardCircle, pos_ba)
    r_b = get_radius(b)
    projection = get_projection(a, pos_ba)
    vec = pos_ba - projection
    normal = LA.normalize(vec)
    if all(isfinite.(normal))
        penetration = r_b - LA.norm(vec)
        contact = pos_ba - (r_b - penetration / 2) * normal
        return Manifold(penetration, normal, contact)
    else
        d, edge_id = get_closest_edge_from_inside(a, pos_ba)
        penetration = r_b + d
        normal = get_normals(a)[edge_id]
        contact = pos_ba - (r_b - penetration / 2) * normal
        return Manifold(penetration, normal, contact)
    end
end

function Manifold(a::StandardCircle, b::StandardRect, pos_ba)
    pos_ab = -pos_ba
    manifold_ab = Manifold(b, a, pos_ab)

    penetration = get_penetration(manifold_ab)
    normal = manifold_ab |> get_normal |> rotate_180
    contact = pos_ba + get_contact(manifold_ab)
    return Manifold(penetration, normal, contact)
end

Manifold(a::StandardRect, b::StandardCircle, pos_ba, dir_ba) = Manifold(a, b, pos_ba)

function Manifold(a::StandardCircle, b::StandardRect, pos_ba, dir_ba)
    dir_ab = invert_relative_direction(dir_ba)
    pos_ab = -rotate(pos_ba, dir_ab)
    manifold_ab = Manifold(b, a, pos_ab, dir_ab)

    penetration = get_penetration(manifold_ab)
    normal = get_relative_direction(dir_ab, rotate_180(get_normal(manifold_ab)))
    contact = rotate(get_contact(manifold_ab) - pos_ab, dir_ba)
    return Manifold(penetration, normal, contact)
end

#####
# StandardRect vs. StandardRect
#####

function Manifold(a::StandardRect, b::StandardRect, pos_ba)
    half_width_a = get_half_width(a)
    half_height_a = get_half_height(a)

    bottom_left = max.(get_bottom_left(a), get_bottom_left(b, pos_ba))
    top_right = min.(get_top_right(a), get_top_right(b, pos_ba))

    contact = (top_right + bottom_left) / 2
    x_contact = contact[1]
    y_contact = contact[2]

    half_widths_intersection = (top_right - bottom_left) / 2
    half_width_intersection = half_widths_intersection[1]
    half_height_intersection = half_widths_intersection[2]

    penetration, edge_id = findmin((half_height_a + y_contact + half_height_intersection,
                                    half_width_a - x_contact + half_width_intersection,
                                    half_height_a - y_contact + half_height_intersection,
                                    half_width_a + x_contact + half_width_intersection,
                                   ))

    normal = get_normals(a)[edge_id]

    return Manifold(penetration, normal, contact)
end

function get_clipped_vertices(a::StandardRect, b::StandardRect, pos_ba, dir_ba)
    half_width_a = get_half_width(a)
    half_height_a = get_half_height(a)

    vertices_ba = get_vertices(b, pos_ba, dir_ba)
    VecType = typeof(pos_ba)

    initial_vertices = (vertices_ba..., vertices_ba[1])
    final_vertices = VecType[]

    # e1
    for i in 1:length(initial_vertices) - 1
        v1 = initial_vertices[i]
        v2 = initial_vertices[i + 1]

        x1 = v1[1]
        y1 = v1[2]

        x2 = v2[1]
        y2 = v2[2]

        # in-in
        if (y1 >= -half_height_a) && (y2 >= -half_height_a)
            push!(final_vertices, v2)
        # in-out
        elseif (y1 >= -half_height_a) && (y2 < -half_height_a)
            x_intersection = x1 + (x2 - x1) * (y1 + half_height_a) / (y1 - y2)
            if isfinite(x_intersection)
                push!(final_vertices, VecType(x_intersection, -half_height_a))
            else
                push!(final_vertices, (v1 + v2) / 2)
            end
        # out-in
        elseif (y1 < -half_height_a) && (y2 >= -half_height_a)
            x_intersection = x1 + (x2 - x1) * (y1 + half_height_a) / (y1 - y2)
            if isfinite(x_intersection)
                push!(final_vertices, VecType(x_intersection, -half_height_a))
                push!(final_vertices, v2)
            else
                push!(final_vertices, (v1 + v2) / 2)
                push!(final_vertices, v2)
            end
        end
    end

    # e2
    initial_vertices = (final_vertices..., final_vertices[1])
    final_vertices = VecType[]

    for i in 1:length(initial_vertices) - 1
        v1 = initial_vertices[i]
        v2 = initial_vertices[i + 1]

        x1 = v1[1]
        y1 = v1[2]

        x2 = v2[1]
        y2 = v2[2]

        # in-in
        if (x1 <= half_width_a) && (x2 <= half_width_a)
            push!(final_vertices, v2)
        # in-out
        elseif (x1 <= half_width_a) && (x2 > half_width_a)
            y_intersection = y1 + (y2 - y1) * (half_width_a - x1) / (x2 - x1)
            if isfinite(y_intersection)
                push!(final_vertices, VecType(half_width_a, y_intersection))
            else
                push!(final_vertices, (v1 + v2) / 2)
            end
        # out-in
        elseif (x1 > half_width_a) && (x2 <= half_width_a)
            y_intersection = y1 + (y2 - y1) * (half_width_a - x1) / (x2 - x1)
            if isfinite(y_intersection)
                push!(final_vertices, VecType(half_width_a, y_intersection))
                push!(final_vertices, v2)
            else
                push!(final_vertices, (v1 + v2) / 2)
                push!(final_vertices, v2)
            end
        end
    end

    initial_vertices = (final_vertices..., final_vertices[1])
    final_vertices = VecType[]

    # e3
    for i in 1:length(initial_vertices) - 1
        v1 = initial_vertices[i]
        v2 = initial_vertices[i + 1]

        x1 = v1[1]
        y1 = v1[2]

        x2 = v2[1]
        y2 = v2[2]

        # in-in
        if (y1 <= half_height_a) && (y2 <= half_height_a)
            push!(final_vertices, v2)
        # in-out
        elseif (y1 <= half_height_a) && (y2 > half_height_a)
            x_intersection = x1 + (x2 - x1) * (half_height_a - y1) / (y2 - y1)
            if isfinite(x_intersection)
                push!(final_vertices, VecType(x_intersection, half_height_a))
            else
                push!(final_vertices, (v1 + v2) / 2)
            end
        # out-in
        elseif (y1 > half_height_a) && (y2 <= half_height_a)
            x_intersection = x1 + (x2 - x1) * (half_height_a - y1) / (y2 - y1)
            if isfinite(x_intersection)
                push!(final_vertices, VecType(x_intersection, half_height_a))
                push!(final_vertices, v2)
            else
                push!(final_vertices, (v1 + v2) / 2)
                push!(final_vertices, v2)
            end
        end
    end

    # e4
    initial_vertices = (final_vertices..., final_vertices[1])
    final_vertices = VecType[]

    for i in 1:length(initial_vertices) - 1
        v1 = initial_vertices[i]
        v2 = initial_vertices[i + 1]

        x1 = v1[1]
        y1 = v1[2]

        x2 = v2[1]
        y2 = v2[2]

        # in-in
        if (x1 >= -half_width_a) && (x2 >= -half_width_a)
            push!(final_vertices, v2)
        # in-out
        elseif (x1 >= -half_width_a) && (x2 < -half_width_a)
            y_intersection = y1 + (y2 - y1) * (x1 + half_width_a) / (x1 - x2)
            if isfinite(y_intersection)
                push!(final_vertices, VecType(-half_width_a, y_intersection))
            else
                push!(final_vertices, (v1 + v2) / 2)
            end
        # out-in
        elseif (x1 < -half_width_a) && (x2 >= -half_width_a)
            y_intersection = y1 + (y2 - y1) * (x1 + half_width_a) / (x1 - x2)
            if isfinite(y_intersection)
                push!(final_vertices, VecType(-half_width_a, y_intersection))
                push!(final_vertices, v2)
            else
                push!(final_vertices, (v1 + v2) / 2)
                push!(final_vertices, v2)
            end
        end
    end

    return final_vertices
end

function get_centroid(vertices::Vararg{Vector2D})
    VecType = typeof(vertices[1])
    T = eltype(vertices[1])
    area = zero(T)
    c_x_num = zero(T)
    c_y_num = zero(T)
    den = zero(T)

    augmented_vertices = (vertices..., vertices[1])

    for i in 1:length(augmented_vertices) - 1
        v_i = augmented_vertices[i]
        v_i_plus_1 = augmented_vertices[i + 1]

        x_i = v_i[1]
        y_i = v_i[2]

        x_i_plus_1 = v_i_plus_1[1]
        y_i_plus_1 = v_i_plus_1[2]

        cross = x_i * y_i_plus_1 - x_i_plus_1 * y_i
        c_x_num += (x_i + x_i_plus_1) * cross
        c_y_num += (y_i + y_i_plus_1) * cross
        den += cross
    end

    c_x = c_x_num / (3 * den)
    c_y = c_y_num / (3 * den)
    if isfinite(c_x) && isfinite(c_y)
        return VecType(c_x, c_y)
    else
        centroid = zero(VecType)
        for vertex in vertices
            centroid = centroid + vertex
        end
        return centroid / length(vertices)
    end
end

function get_contact(a::StandardRect, b::StandardRect, pos_ba, dir_ba)
    clipped_vertices = get_clipped_vertices(a, b, pos_ba, dir_ba)
    return get_centroid(clipped_vertices...)
end

function get_candidate_support(a::StandardRect, b::StandardRect, pos_ba, dir_ba)
    half_width_a = get_half_width(a)
    half_height_a = get_half_height(a)

    vertices_ba = get_vertices(b, pos_ba, dir_ba)

    value_1, _, vertex_id_1 = find_maximum(vertex -> vertex[2], vertices_ba)
    max_penetration_1 = half_height_a + value_1

    value_2, _, vertex_id_2 = find_minimum(vertex -> vertex[1], vertices_ba)
    max_penetration_2 = half_width_a - value_2

    value_3, _, vertex_id_3 = find_minimum(vertex -> vertex[2], vertices_ba)
    max_penetration_3 = half_height_a - value_3

    value_4, _, vertex_id_4 = find_maximum(vertex -> vertex[1], vertices_ba)
    max_penetration_4 = half_width_a + value_4

    penetration, (_, vertex_id), edge_id = find_minimum(x -> x[1], ((max_penetration_1, vertex_id_1),
                                          (max_penetration_2, vertex_id_2),
                                          (max_penetration_3, vertex_id_3),
                                          (max_penetration_4, vertex_id_4)))

    return penetration, vertex_id, edge_id
end

function Manifold(a::StandardRect, b::StandardRect, pos_ba, dir_ba)
    penetration_ba, vertex_id_b, edge_id_a = get_candidate_support(a, b, pos_ba, dir_ba)

    pos_ab, dir_ab = invert(pos_ba, dir_ba)
    penetration_ab, vertex_id_a, edge_id_b = get_candidate_support(b, a, pos_ab, dir_ab)

    contact = get_contact(a, b, pos_ba, dir_ba)

    if penetration_ba <= penetration_ab
        penetration = penetration_ba
        normal = get_normals(a)[edge_id_a]
        return Manifold(penetration, normal, contact)
    else
        penetration = penetration_ab
        normal = get_normals(b, dir_ba)[edge_id_b]
        return Manifold(penetration, normal, contact)
    end
end
