struct Manifold{T}
    penetration::T
    axes::Axes{T}
    contact::GB.Vec{2, T}
end

get_penetration(manifold::Manifold) = manifold.penetration
get_axes(manifold::Manifold) = manifold.axes
get_normal(manifold::Manifold) = manifold |> get_axes |> get_y_cap
get_tangent(manifold::Manifold) = manifold |> get_axes |> get_x_cap
get_contact(manifold::Manifold) = manifold.contact

#####
# StdCircle vs. StdCircle
#####

function Manifold(a::StdCircle{T}, b::StdCircle{T}, pos_ba::GB.Vec{2, T}) where {T}
    r_a = get_radius(a)
    r_b = get_radius(b)
    normal = LA.normalize(pos_ba)
    if all(isfinite.(normal))
        penetration = r_a + r_b - LA.norm(pos_ba)
        tangent = rotate_minus_90(normal)
        axes = Axes(tangent, normal)
        contact = (r_a - penetration / 2) .* normal
        return Manifold(penetration, axes, contact)
    else
        penetration = r_a + r_b
        axes = rotate_minus_90(Axes{T}())
        normal = get_x_cap(axes)
        contact = (r_a - penetration / 2) .* normal
        return Manifold(penetration, axes, contact)
    end
end

Manifold(a::StdCircle{T}, b::StdCircle{T}, pos_ba::GB.Vec{2, T}, axes_ba::Axes{T}) where {T} = Manifold(a, b, pos_ba)

#####
# HyperRectangle vs. HyperSphere
#####

function get_closest_edge_from_inside(rect::StdRect{T}, pos::GB.Vec{2, T}) where {T}
    half_width = get_half_width(rect)
    half_height = get_half_height(rect)

    x = get_x(pos)
    y = get_y(pos)

    penetration, edge_id = findmin((y + half_height, half_width - x, half_height - y, x + half_width))

    return (penetration, edge_id)
end

function Manifold(a::StdRect{T}, b::StdCircle{T}, pos_ba::GB.Vec{2, T}) where {T}
    r_b = get_radius(b)
    projection = get_projection(a, pos_ba)
    vec = pos_ba .- projection
    normal = LA.normalize(vec)
    if all(isfinite.(normal))
        penetration = r_b - LA.norm(vec)
        tangent = rotate_minus_90(normal)
        axes = Axes(tangent, normal)
        contact = pos_ba .+ (r_b - penetration / 2) .* -normal
        return Manifold(penetration, axes, contact)
    else
        d, edge_id = get_closest_edge_from_inside(a, pos_ba)
        penetration = r_b + d
        normal = get_normals(a)[edge_id]
        tangent = rotate_minus_90(normal)
        axes = Axes(tangent, normal)
        contact = pos_ba .+ (r_b - penetration / 2) .* -normal
        return Manifold(penetration, axes, contact)
    end
end

function Manifold(a::StdCircle{T}, b::StdRect{T}, pos_ba::GB.Vec{2, T}) where {T}
    pos_ab = -pos_ba
    manifold_ab = Manifold(b, a, pos_ab)

    penetration = get_penetration(manifold_ab)
    axes = manifold_ab |> get_axes |> rotate_180
    contact = pos_ba .+ get_contact(manifold_ab)
    return Manifold(penetration, axes, contact)
end

Manifold(a::StdRect{T}, b::StdCircle{T}, pos_ba::GB.Vec{2, T}, axes_ba::Axes{T}) where {T} = Manifold(a, b, pos_ba)

function Manifold(a::StdCircle{T}, b::StdRect{T}, pos_ba::GB.Vec{2, T}, axes_ba::Axes{T}) where {T}
    axes_ab = invert_relative_axes(axes_ba)
    pos_ab = -rotate(pos_ba, axes_ab)
    manifold_ab = Manifold(b, a, pos_ab, axes_ab)

    penetration = get_penetration(manifold_ab)
    axes = get_relative_axes(axes_ab, rotate_180(get_axes(manifold_ab)))
    contact = rotate(get_contact(manifold_ab) .- pos_ab, axes_ba)
    return Manifold(penetration, axes, contact)
end

#####
# StdRect vs. StdRect
#####

function Manifold(a::StdRect{T}, b::StdRect{T}, pos_ba::GB.Vec{2, T}) where {T}
    half_width_a = get_half_width(a)
    half_height_a = get_half_height(a)

    bottom_left = max.(get_bottom_left(a), get_bottom_left(b, pos_ba))
    top_right = min.(get_top_right(a), get_top_right(b, pos_ba))

    contact = (top_right .+ bottom_left) ./ 2
    x_contact = get_x(contact)
    y_contact = get_y(contact)

    half_widths_intersection = (top_right .- bottom_left) ./ 2
    half_width_intersection = half_widths_intersection[1]
    half_height_intersection = half_widths_intersection[2]

    penetration, edge_id = findmin((half_height_a + y_contact + half_height_intersection,
                                    half_width_a - x_contact + half_width_intersection,
                                    half_height_a - y_contact + half_height_intersection,
                                    half_width_a + x_contact + half_width_intersection,
                                   ))

    normal = get_normals(a)[edge_id]
    tangent = rotate_minus_90(normal)
    axes = Axes(tangent, normal)

    return Manifold(penetration, axes, contact)
end

function get_clipped_vertices(a::StdRect{T}, b::StdRect{T}, pos_ba::GB.Vec{2, T}, axes_ba::Axes{T}) where {T}
    half_width_a = get_half_width(a)
    half_height_a = get_half_height(a)

    vertices_ba = get_vertices(b, pos_ba, axes_ba)
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

function get_centroid(vertices::Vararg{GB.Vec{2, T}}) where {T}
    VecType = typeof(vertices[1])
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
            centroid = centroid .+ vertex
        end
        return centroid / length(vertices)
    end
end

function get_contact(a::StdRect{T}, b::StdRect{T}, pos_ba::GB.Vec{2, T}, axes_ba::Axes{T}) where {T}
    clipped_vertices = get_clipped_vertices(a, b, pos_ba, axes_ba)
    return get_centroid(clipped_vertices...)
end

function get_candidate_support(a::StdRect{T}, b::StdRect{T}, pos_ba::GB.Vec{2, T}, axes_ba::Axes{T}) where {T}
    half_width_a = get_half_width(a)
    half_height_a = get_half_height(a)

    vertices_ba = get_vertices(b, pos_ba, axes_ba)

    value_1, _, vertex_id_1 = findmax(vertex -> vertex[2], vertices_ba)
    max_penetration_1 = half_height_a + value_1

    value_2, _, vertex_id_2 = findmin(vertex -> vertex[1], vertices_ba)
    max_penetration_2 = half_width_a - value_2

    value_3, _, vertex_id_3 = findmin(vertex -> vertex[2], vertices_ba)
    max_penetration_3 = half_height_a - value_3

    value_4, _, vertex_id_4 = findmax(vertex -> vertex[1], vertices_ba)
    max_penetration_4 = half_width_a + value_4

    penetration, (_, vertex_id), edge_id = findmin(x -> x[1], ((max_penetration_1, vertex_id_1),
                                          (max_penetration_2, vertex_id_2),
                                          (max_penetration_3, vertex_id_3),
                                          (max_penetration_4, vertex_id_4)))

    return penetration, vertex_id, edge_id
end

function Manifold(a::StdRect{T}, b::StdRect{T}, pos_ba::GB.Vec{2, T}, axes_ba::Axes{T}) where {T}
    penetration_ba, vertex_id_b, edge_id_a = get_candidate_support(a, b, pos_ba, axes_ba)

    pos_ab, axes_ab = invert(pos_ba, axes_ba)
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
