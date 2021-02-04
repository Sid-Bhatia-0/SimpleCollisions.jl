const ATOL = 1e-7

function is_colliding(a::RigidBody, b::RigidBody)
    shape_a = get_shape(a)
    shape_b = get_shape(b)

    pos_a = get_position(a)
    pos_b = get_position(b)
    pos_ba = position_b .- position_a

    axes_a = get_axes(a)
    axes_b = get_axes(b)
    axes_ba = get_relative_axes(axes_a, axes_b)

    is_colliding(shape_a, shape_b, pos_ba, axes_ba)
end

#####
# Point vs. Point
#####

is_colliding(a::GB.Vec, b::GB.Vec, pos_ba) = pos_ba == zero(pos_ba)
is_colliding(a::GB.Vec, b::GB.Vec, pos_ba, axes_ba) = is_colliding(a, b, pos_ba)

#####
# Line vs. Point
#####

function is_colliding(a::GB.Line, b::GB.Vec, pos_ba)
    p1 = get_point(a, 1)
    p2 = get_point(a, 2)
    x_ba = pos_ba[1]
    y_ba = pos_ba[2]
    return (y_ba == zero(y_ba)) && (p1[1] <= x_ba <= p2[1])
end

is_colliding(a::GB.Line, b::GB.Vec, pos_ba, axes_ba) = is_colliding(a, b, pos_ba)
function is_colliding(a::GB.Vec, b::GB.Line, pos_ba, axes_ba)
    axes_ab = invert_relative_axes(axes_ba)
    pos_ab = -rotate(pos_ba, axes_ab)
    is_colliding(b, a, pos_ab, axes_ab)
end

#####
# Line vs. Line
#####

function is_colliding(a::GB.Line, b::GB.Line, pos_ba, axes_ba)
    p1 = get_point(a, 1)
    p2 = get_point(a, 2)

    half_width_b = get_point(b, 2)[1]
    x_cap_b = get_x_cap(axes_ba)
    q1 = pos_ba .- half_width_b .* x_cap_b
    q2 = pos_ba .+ half_width_b .* x_cap_b

    x1_a = p1[1]
    x2_a = p2[1]

    x1_b = q1[1]
    y1_b = q1[2]
    x2_b = q2[1]
    y2_b = q2[2]

    t = y1_b / (y1_b - y2_b)

    if isfinite(t)
        x_intersection = x1_b + t * (x2_b - x1_b)
        return (zero(t) <= t <= one(t)) && (x1_a <= x_intersection <= x2_a)
    else
        return (y1_b == zero(y1_b)) && (y2_b == zero(y2_b)) && !((x2_a < x1_b) || (x2_b < x1_a))
    end
end

#####
# HyperSphere vs. Point
#####

is_colliding(a::GB.HyperSphere, b::GB.Vec, pos_ba) = LA.dot(pos_ba, pos_ba) <= a.r ^ 2
is_colliding(a::GB.HyperSphere, b::GB.Vec, pos_ba, axes_ba) = is_colliding(a, b, pos_ba)
is_colliding(a::GB.Vec, b::GB.HyperSphere, pos_ba, axes_ba) = is_colliding(b, a, pos_ba)

#####
# HyperSphere vs. Line
#####

function is_colliding(a::GB.HyperSphere, b::GB.Line, pos_ba, axes_ba)
    axes_ab = invert_relative_axes(axes_ba)
    pos_ab = -rotate(pos_ba, axes_ab)
    is_colliding(b, a, pos_ab, axes_ab)
end

is_colliding(a::GB.Line, b::GB.HyperSphere, pos_ba, axes_ba) = is_colliding(a, b, pos_ba)

function project(a::GB.Line, b::GB.HyperSphere, pos_ba)
    VecType = typeof(pos_ba)
    half_width_a = get_point(a, 2)[1]
    x_b = pos_ba[1]
    if x_b < -half_width_a
        return get_point(a, 1)
    elseif x_b > half_width_a
        return get_point(a, 2)
    else
        return VecType(x_b, zero(x_b))
    end
end

function is_colliding(a::GB.Line, b::GB.HyperSphere, pos_ba)
    closest_point_ba = project(a, b, pos_ba)
    vec = pos_ba .- closest_point_ba
    return LA.dot(vec, vec) <= b.r ^ 2
end

#####
# HyperSphere vs. HyperSphere
#####

is_colliding(a::GB.HyperSphere, b::GB.HyperSphere, pos_ba) = LA.dot(pos_ba, pos_ba) <= (a.r + b.r) ^ 2
is_colliding(a::GB.HyperSphere, b::GB.HyperSphere, pos_ba, axes_ba) = is_colliding(a, b, pos_ba)

#####
# HyperRectangle vs. Point
#####

is_colliding(a::GB.Rect, b::GB.Vec, pos_ba) = all(minimum(a) .<= pos_ba .<= maximum(a))
is_colliding(a::GB.Rect, b::GB.Vec, pos_ba, axes_ba) = is_colliding(a, b, pos_ba)
function is_colliding(a::GB.Vec, b::GB.Rect, pos_ba, axes_ba)
    axes_ab = invert_relative_axes(axes_ba)
    pos_ab = -rotate(pos_ba, axes_ab)
    return is_colliding(b, a, pos_ab, axes_ab)
end

#####
# HyperRectangle vs. Line
#####

function is_colliding(a::GB.Rect, b::GB.Line, pos_ba, axes_ba)
    origin = zero(a.origin)
    half_width_b = get_point(b, 2)[1]
    x_cap_ba = get_x_cap(axes_ba)
    q1 = pos_ba .+ half_width_b .* x_cap_ba
    q2 = pos_ba .- half_width_b .* x_cap_ba

    if is_colliding(a, origin, q1) || is_colliding(a, origin, q2)
        return true
    else
        half_widths_a = get_half_widths(a)
        half_width_a = half_widths_a[1]
        half_height_a = half_widths_a[2]

        PointType = typeof(b.points[1])
        VecType = typeof(pos_ba)

        x_ba = pos_ba[1]
        y_ba = pos_ba[2]

        zero_val = zero(half_width_a)
        horizontal_line = GB.Line(PointType(-half_width_a, zero_val), PointType(half_width_a, zero_val))
        vertical_line = GB.Line(PointType(-half_height_a, zero_val), PointType(half_height_a, zero_val))

        axis_ba_minus_90 = rotate_minus_90(axes_ba)
        return any([is_colliding(horizontal_line, b, VecType(x_ba, y_ba + half_height_a), axes_ba),
                    is_colliding(vertical_line, b, VecType(y_ba, half_width_a - x_ba), axis_ba_minus_90),
                    is_colliding(horizontal_line, b, VecType(x_ba, y_ba - half_height_a), axes_ba),
                    is_colliding(vertical_line, b, VecType(y_ba, - half_width_a - x_ba), axis_ba_minus_90)])
    end
end

function is_colliding(a::GB.Line, b::GB.Rect, pos_ba, axes_ba)
    axes_ab = invert_relative_axes(axes_ba)
    pos_ab = -rotate(pos_ba, axes_ab)
    return is_colliding(b, a, pos_ab, axes_ab)
end

#####
# HyperRectangle vs. HyperSphere
#####

project(a::GB.Rect, b::GB.HyperSphere, pos_ba) = clamp.(pos_ba, minimum(a), maximum(a))

function is_colliding(a::GB.Rect, b::GB.HyperSphere, pos_ba)
    closest_point_ba = project(a, b, pos_ba)
    vec = pos_ba .- closest_point_ba
    return LA.dot(vec, vec) <= b.r ^ 2
end

is_colliding(a::GB.Rect, b::GB.HyperSphere, pos_ba, axes_ba) = is_colliding(a, b, pos_ba)
is_colliding(a::GB.HyperSphere, b::GB.Rect, pos_ba, axes_ba) = is_colliding(b, a, -pos_ba)

#####
# HyperRectangle vs. HyperRectangle
#####

function is_penetrating(a::GB.Rect, b::GB.Rect, pos_ba, axes_ba)
    half_widths_a = maximum(a)
    half_width_a = half_widths_a[1]
    half_height_a = half_widths_a[2]

    bottom_left_ba, bottom_right_ba, top_right_ba, top_left_ba = get_vertices(b, pos_ba, axes_ba)

    min_x_ba, max_x_ba = extrema((bottom_left_ba[1], bottom_right_ba[1], top_right_ba[1], top_left_ba[1]))
    min_y_ba, max_y_ba = extrema((bottom_left_ba[2], bottom_right_ba[2], top_right_ba[2], top_left_ba[2]))

    return !((half_width_a <= min_x_ba) || (max_x_ba <= -half_width_a) || (half_height_a <= min_y_ba) || (max_y_ba <= -half_height_a))
end

is_colliding(a::GB.Rect, b::GB.Rect, pos_ba, axes_ba) = is_penetrating(a, b, pos_ba, axes_ba) || is_penetrating(b, a, -pos_ba, invert_relative_axes(axes_ba))
