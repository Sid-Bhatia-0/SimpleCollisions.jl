struct Manifold{T}
    penetration::T
    axes::Axes{T}
    contact::GB.Vec2{T}
end

get_penetration(manifold::Manifold) = manifold.penetration
get_axes(manifold::Manifold) = manifold.axes
get_normal(manifold::Manifold) = manifold |> get_axes |> get_y_cap
get_tangent(manifold::Manifold) = manifold |> get_axes |> get_x_cap
get_contact(manifold::Manifold) = manifold.contact

function Manifold(a::RigidBody, b::RigidBody)
    shape_a = get_shape(a)
    shape_b = get_shape(b)

    pos_a = get_position(a)
    pos_b = get_position(b)
    pos_ba = position_b .- position_a

    axes_a = get_axes(a)
    axes_b = get_axes(b)
    axes_ba = get_relative_axes(axes_a, axes_b)

    return Manifold(shape_a, shape_b, pos_ba, axes_ba)
end

#####
# HyperSphere vs. HyperSphere
#####

function Manifold(a::GB.HyperSphere, b::GB.HyperSphere, pos_ba)
    normal = LA.normalize(pos_ba)
    if all(isfinite.(normal))
        d = LA.norm(pos_ba)
        penetration = a.r + b.r - d
        tangent = rotate_minus_90(normal)
        axes = Axes(tangent, normal)
        contact = (a.r - penetration / 2) .* normal
        return Manifold(penetration, axes, contact)
    else
        penetration = a.r + b.r
        T = eltype(pos_ba)
        axes = rotate_minus_90(Axes{T}())
        tangent = get_x_cap(axes)
        normal = get_y_cap(axes)
        contact = (a.r - penetration / 2) .* normal
        return Manifold(penetration, axes, contact)
    end
end

Manifold(a::GB.HyperSphere, b::GB.HyperSphere, pos_ba, axes_ba) = Manifold(a, b, pos_ba)

#####
# HyperRectangle vs. HyperSphere
#####

function project_from_inside(a::GB.HyperRectangle{N}, b::GB.Vec{N}) where {N}
    bottom_left = get_bottom_left(a)
    top_right = get_top_right(a)
    min_dist = Inf
    min_dist_dim = 0
    min_dist_dir = 0

    for i in 1:N
        dist, dir = findmin((b[i] - bottom_left[i], top_right[i] - b[i]))
        if dist < min_dist
            min_dist = dist
            min_dist_dim = i
            min_dist_dir = dir == 1 ? -1 : 1
        end
    end

    return (min_dist, min_dist_dim, min_dist_dir)
end

function Manifold(a::GB.Rect, b::GB.HyperSphere, pos_ba::GB.Vec)
    closest_point_ba = project(a, b, pos_ba)
    vec = pos_ba .- closest_point_ba
    normal = LA.normalize(vec)
    if all(isfinite.(normal))
        d = LA.norm(vec)
        penetration = b.r - d
        tangent = rotate_minus_90(normal)
        axes = Axes(tangent, normal)
        contact = pos_ba .+ (b.r - penetration / 2) .* -normal
        return Manifold(penetration, axes, contact)
    else
        dist, dim, dir = project_from_inside(a, pos_ba)
        penetration = b.r + dist
        normal = dir * GB.unit(typeof(pos_ba), dim)
        tangent = rotate_minus_90(normal)
        axes = Axes(tangent, normal)
        contact = pos_ba .+ (b.r - penetration / 2) .* -normal
        return Manifold(penetration, axes, contact)
    end
end

Manifold(a::GB.Rect, b::GB.HyperSphere, pos_ba, axes_ba) = Manifold(a, b, pos_ba)

function Manifold(a::GB.HyperSphere, b::GB.Rect, pos_ba, axes_ba)
    axes_ab = invert_relative_axes(axes_ba)
    pos_ab = -rotate(pos_ba, axes_ab)
    manifold_ab = Manifold(b, a, pos_ab, axes_ab)

    penetration = get_penetration(manifold_ab)
    axes = get_relative_axes(axes_ab, rotate_180(get_axes(manifold_ab)))
    contact = rotate(get_contact(manifold_ab) .- pos_ab, axes_ba)
    return Manifold(penetration, axes, contact)
end

#####
# HyperRectangle vs. HyperRectangle
#####

function get_clipped_vertices(a::GB.Rect, b::GB.Rect, pos_ba, axes_ba)
    half_widths_a = get_top_right(a)
    half_width_a = half_widths_a[1]
    half_height_a = half_widths_a[2]

    vertices_ba = get_vertices(b, pos_ba, axes_ba)
    VecType = typeof(vertices_ba[1])

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
                push!(final_vertices, (v1 .+ v2) ./ 2)
            end
        # out-in
        elseif (y1 < -half_height_a) && (y2 >= -half_height_a)
            x_intersection = x1 + (x2 - x1) * (y1 + half_height_a) / (y1 - y2)
            if isfinite(x_intersection)
                push!(final_vertices, VecType(x_intersection, -half_height_a))
                push!(final_vertices, v2)
            else
                push!(final_vertices, (v1 .+ v2) ./ 2)
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
                push!(final_vertices, (v1 .+ v2) ./ 2)
            end
        # out-in
        elseif (x1 > half_width_a) && (x2 <= half_width_a)
            y_intersection = y1 + (y2 - y1) * (half_width_a - x1) / (x2 - x1)
            if isfinite(y_intersection)
                push!(final_vertices, VecType(half_width_a, y_intersection))
                push!(final_vertices, v2)
            else
                push!(final_vertices, (v1 .+ v2) ./ 2)
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
                push!(final_vertices, (v1 .+ v2) ./ 2)
            end
        # out-in
        elseif (y1 > half_height_a) && (y2 <= half_height_a)
            x_intersection = x1 + (x2 - x1) * (half_height_a - y1) / (y2 - y1)
            if isfinite(x_intersection)
                push!(final_vertices, VecType(x_intersection, half_height_a))
                push!(final_vertices, v2)
            else
                push!(final_vertices, (v1 .+ v2) ./ 2)
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
                push!(final_vertices, (v1 .+ v2) ./ 2)
            end
        # out-in
        elseif (x1 < -half_width_a) && (x2 >= -half_width_a)
            y_intersection = y1 + (y2 - y1) * (x1 + half_width_a) / (x1 - x2)
            if isfinite(y_intersection)
                push!(final_vertices, VecType(-half_width_a, y_intersection))
                push!(final_vertices, v2)
            else
                push!(final_vertices, (v1 .+ v2) ./ 2)
                push!(final_vertices, v2)
            end
        end
    end

    return final_vertices
end

function get_centroid(augmented_vertices...)
    VecType = typeof(augmented_vertices[1])
    T = eltype(augmented_vertices[1])
    area = zero(T)
    c_x_num = zero(T)
    c_y_num = zero(T)
    den = zero(T)

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
        return reduce(+, augmented_vertices[1:end-1]) ./ (length(augmented_vertices) - 1)
    end
end

function get_contact(a::GB.Rect, b::GB.Rect, pos_ba, axes_ba)
    clipped_vertices = get_clipped_vertices(a, b, pos_ba, axes_ba)
    return get_centroid(clipped_vertices..., clipped_vertices[1])
end

function get_candidate_support(a::GB.Rect, b::GB.Rect, pos_ba, axes_ba)
    half_widths_a = get_top_right(a)
    half_width_a = half_widths_a[1]
    half_height_a = half_widths_a[2]

    vertices_ba = get_vertices(b, pos_ba, axes_ba)

    value_1, vertex_1, vertex_id_1 = findmax(vertex -> vertex[2], vertices_ba)
    max_penetration_1 = half_height_a + value_1

    value_2, vertex_2, vertex_id_2 = findmin(vertex -> vertex[1], vertices_ba)
    max_penetration_2 = half_width_a - value_2

    value_3, vertex_3, vertex_id_3 = findmin(vertex -> vertex[2], vertices_ba)
    max_penetration_3 = half_height_a - value_3

    value_4, vertex_4, vertex_id_4 = findmax(vertex -> vertex[1], vertices_ba)
    max_penetration_4 = half_width_a + value_4

    penetration, (max_penetration, vertex_id), edge_id = findmin(x -> x[1], ((max_penetration_1, vertex_id_1),
                                          (max_penetration_2, vertex_id_2),
                                          (max_penetration_3, vertex_id_3),
                                          (max_penetration_4, vertex_id_4)))

    return penetration, vertex_id, edge_id
end

function Manifold(a::GB.Rect, b::GB.Rect, pos_ba, axes_ba)
    penetration_ba, vertex_id_b, edge_id_a = get_candidate_support(a, b, pos_ba, axes_ba)

    axes_ab = invert_relative_axes(axes_ba)
    pos_ab = -rotate(pos_ba, axes_ab)
    penetration_ab, vertex_id_a, edge_id_b = get_candidate_support(b, a, pos_ab, axes_ab)

    contact = get_contact(a, b, pos_ba, axes_ba)

    if penetration_ba <= penetration_ab
        penetration = penetration_ba
        normal = get_normals(a)[edge_id_a]
        tangent = rotate_minus_90(normal)
        axes = Axes(tangent, normal)
        return Manifold(penetration, axes, contact)
    else
        penetration = penetration_ab
        normal = get_normals(b, pos_ba, axes_ba)[edge_id_b]
        tangent = rotate_minus_90(normal)
        axes = Axes(tangent, normal)
        return Manifold(penetration, axes, contact)
    end
end
