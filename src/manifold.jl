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

function get_relative_manifold(manifold_ab::Manifold, pos_ba, axes_ba)
    penetration = get_penetration(manifold_ab)
    axes = get_relative_axes(axes_ba, rotate_180(get_axes(manifold_ab)))
    contact = get_contact(manifold_ab) .- pos_ba
    return Manifold(penetration, axes, contact)
end

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
        axes = Axes{T}()
        contact = zero(pos_ba)
        return Manifold(penetration, normal, contact)
    end
end

#####
# HyperRectangle vs. HyperSphere
#####

function project_edge(a::GB.HyperRectangle{N}, b::GB.Vec{N}) where {N}
    bottom_left = minimum(a)
    top_right = maximum(a)
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

function Manifold(a::GB.Rect, b::GB.HyperSphere, pos_ba)
    closest_point_ba = project(a, b, pos_ba)
    vec = pos_ba .- closest_point_ba
    normal = LA.normalize(vec)
    if all(isfinite.(normal))
        d = LA.norm(vec)
        penetration = b.r - d
        tangent = rotate_minus_90(normal)
        axes = Axes(tangent, normal)
        contact = pos_ba .+ (a.r - penetration / 2) .* -normal
        return Manifold(penetration, axes, contact)
    else
        dist, dim, dir = project_edge(a, pos_ba)
        penetration = b.r + dist
        normal = dir * GB.unit(typeof(pos_ba), dim)
        tangent = rotate_minus_90(normal)
        axes = Axes(tangent, normal)
        contact = pos_ba .+ (a.r - penetration / 2) .* -normal
        return Manifold(penetration, axes, contact)
    end
end

function Manifold(a::GB.HyperSphere, b::GB.Rect, pos_ba, axes_ba)
    manifold_ab = Manifold(b, a, -pos_ba)
    manifold_ba = get_relative_manifold(manifold_ab, pos_ba, axes_ba)
    return manifold_ba
end

#####
# HyperRectangle vs. HyperRectangle
#####

function get_candidate_support(a::GB.Rect, b::GB.Rect, pos_ba, axes_ba)
    half_widths_a = maximum(a)
    half_width_a = half_widths_a[1]
    half_height_a = half_widths_a[2]

    vertices_ba = get_vertices(b, pos_ba, axes_ba)

    max_value_1, max_vertex_1, max_vertex_id_1 = findmax(vertex -> vertex[2], vertices_ba)
    max_penetration_1 = half_height_a + max_value_1

    min_value_2, max_vertex_2, max_vertex_id_2 = findmin(vertex -> vertex[1], vertices_ba)
    max_penetration_2 = half_width_a - min_value_2

    min_value_3, max_vertex_3, max_vertex_id_3 = findmin(vertex -> vertex[2], vertices_ba)
    max_penetration_3 = half_height_a - min_value_3

    max_value_4, max_vertex_4, max_vertex_id_4 = findmax(vertex -> vertex[1], vertices_ba)
    max_penetration_4 = half_width_a + max_value_4

    max_penetration, candidate_support, edge_id = findmin(x -> x[1], ((max_penetration_1, max_vertex_1, max_vertex_id_1),
                                          (max_penetration_2, max_vertex_2, max_vertex_id_2),
                                          (max_penetration_3, max_vertex_3, max_vertex_id_3),
                                          (max_penetration_4, max_vertex_4, max_vertex_id_4)))

    max_vertex = candidate_support[2]
    vertex_id = candidate_support[3]

    return max_penetration, vertex, vertex_id, edge_id
end

function Manifold(a::GB.Rect, b::GB.Rect, pos_ba, axes_ba)
    penetration_ba, vertex_ba, vertex_id_ba, edge_id_ba = get_candidate_support(a, b, pos_ba, axes_ba)
    penetration_ab, vertex_ab, vertex_id_ab, edge_id_ab = get_candidate_support(b, a, -pos_ba, invert_relative_axes(axes_ba))
    if penetration_ba <= penetration_ab
        penetration = penetration_ba
        normals_aa = get_normals(a)
        normal = normals_aa[edge_id_ba]
        tangent = rotate_minus_90(normal)
        axes = Axes(tangent, normal)
        contact = vertex_ba .+ (penetration / 2) .* -normal
        return Manifold(penetration, axes, contact)
    else
        penetration = penetration_ab
        normals_ba = get_normals(b, pos_ba, axes_ba)
        normal = normals_ba[edge_id_ab]
        tangent = rotate_minus_90(normal)
        axes = Axes(tangent, normal)
        vertices_ba = get_vertices(b, pos_ba, axes_ba)
        vertex = vertices_ba[vertex_id_ab]
        contact = vertex .+ (penetration / 2) .* -normal
        return Manifold(penetration, axes, contact)
    end
end
