const ATOL = 1e-7

get_relative_direction(d1::GB.Vec, d2::GB.Vec) = typeof(d1)(d2[1] * d1[1] + d2[2] * d1[2], d2[2] * d1[1] - d2[1] * d1[2])
get_reverse_direction(d::GB.Vec) = typeof(d)(d[1], -d[2])
get_minus_90_rotated(d::GB.Vec) = typeof(d)(d[2], -d[1])
get_plus_90_rotated(d::GB.Vec) = typeof(d)(-d[2], d[1])

function is_colliding(a::RigidBody, b::RigidBody)
    shape_a = get_shape(a)
    shape_b = get_shape(b)

    position_a = get_position(a)
    position_b = get_position(b)
    position_ba = position_b .- position_a

    direction_a = get_direction(a)
    direction_b = get_direction(b)
    direction_ba = get_relative_direction(direction_a, direction_b)

    is_colliding(shape_a, shape_b, position_ba, direction_ba)
end

#####
# Point vs. Point
#####

is_colliding(a::GB.Vec, b::GB.Vec, position_ba, direction_ba) = all(ab -> isapprox(ab[1], ab[2], atol = ATOL), zip(zero(position_ba), position_ba))
# is_colliding(a::GB.Vec{N}, b::GB.Vec{N}) where {N} = all(ab -> isapprox(ab[1], ab[2], atol = ATOL), zip(a, b))

#####
# Line vs. Point
#####

function is_colliding(a::GB.Line, b::GB.Vec, position_ba, direction_ba)
    p1 = get_point(a, 1)
    p2 = get_point(a, 2)
    x = position_ba[1]
    y = position_ba[2]
    return isapprox(y, zero(y), atol = ATOL) && (p1[1] <= x <= p2[1])
end

is_colliding(a::GB.Vec, b::GB.Line, position_ba, direction_ba) = is_colliding(b, a, -position_ba, get_reverse_direction(direction_ba))

# function is_colliding(a::GB.Line{N}, b::GB.Vec{N}) where {N}
    # T = eltype(b)
    # p1 = a.points[1]
    # p2 = a.points[2]
    # p1_b = p1 .- b
    # p2_b = p2 .- b
    # return LA.dot(p1_b, p2_b) <= zero(T) && isapprox(get_area(b, p1, p2), zero(T), atol = ATOL)
# end

# is_colliding(a::GB.Vec{N}, b::GB.Line{N}) where {N} = is_colliding(b, a)

#####
# Line vs. Line
#####

function is_colliding(a::GB.Line, b::GB.Line, position_ba, direction_ba)
    p1 = get_point(a, 1)
    p2 = get_point(a, 2)
    norm_b = get_point(b, 1)[1] - get_point(b, 2)[1]
    q1 = position_ba .+ (norm_b / 2) * direction_ba
    q2 = position_ba .- (norm_b / 2) * direction_ba

    x1 = q1[1]
    y1 = q1[2]
    x2 = q2[1]
    y2 = q2[2]
    if y1 <= zero(y1) <= y2
        if isapprox(y1, y2, atol = ATOL)
            return p1[1] <= x1 <= p2[1] || p1[1] <= x2 <= p2[1]
        else
            x1 = q1[1]
            x2 = q2[1]
            x = x1 + (y1 / (y1 - y2)) * (x2 - x1)
            return p1[1] <= x <= p2[1]
        end
    else
        return false
    end
end

# function is_colliding(a::GB.Line{2}, b::GB.Line{2})
    # p1 = get_point(a, 1)
    # p2 = get_point(a, 2)
    # q1 = get_point(b, 1)
    # q2 = get_point(b, 2)
    # p2_p1 = p2 .- p1
    # q2_q1 = q2 .- q1

    # normal_a = GB.Vec(-p2_p1[2], p2_p1[1])
    # x = LA.dot(q2_q1, normal_a)
    # if isapprox(x, zero(x), atol = ATOL)
        # return is_colliding(a, q1) || is_colliding(a, q2)
    # else
        # normal_b = GB.Vec(-q2_q1[2], q2_q1[1])
        # y = LA.dot(p2_p1, normal_b)
        # g = LA.dot(q1 .- p1, normal_b) / y # note that since the lines are not parallel (since x is non-zero), y will not be zero
        # h = LA.dot(p1 .- q1, normal_a) / x
        # if zero(g) <= g <= one(g) && zero(h) <= h <= one(h)
            # return true
        # else
            # return false
        # end
    # end
# end

#####
# HyperSphere vs. Point
#####

is_colliding(a::GB.HyperSphere, b::GB.Vec, position_ba, direction_ba) = LA.dot(position_ba, position_ba) <= a.r ^ 2

is_colliding(a::GB.Vec, b::GB.HyperSphere, position_ba, direction_ba) = is_colliding(b, a, -position_ba, get_reverse_direction(direction_ba))

# function is_colliding(a::GB.HyperSphere{N}, b::GB.Vec{N}) where {N}
    # center_a = get_center(a)
    # ba = b .- center_a
    # return LA.dot(ba, ba) <= a.r ^ 2
# end

# is_colliding(a::GB.Vec{N}, b::GB.HyperSphere{N}) where {N} = is_colliding(b, a)

#####
# HyperSphere vs. Line
#####

is_colliding(a::GB.HyperSphere, b::GB.Line, position_ba, direction_ba) = is_colliding(b, a, -position_ba, get_reverse_direction(direction_ba))

function is_colliding(a::GB.Line, b::GB.HyperSphere, position_ba, direction_ba)
end

# function is_colliding(a::GB.HyperSphere{N}, b::GB.Line{N}) where {N}
    # p1 = get_point(b, 1)
    # p2 = get_point(b, 2)
    # if is_colliding(a, p1) || is_colliding(a, p2)
        # return true
    # else
        # center_a = get_center(a)
        # area = get_area(center_a, p1, p2)
        # p2_p1 = p2 .- p1
        # height_squared = 4 * area ^ 2 / LA.dot(p2_p1, p2_p1)
        # if height_squared > a.r ^ 2
            # return false
        # else
            # a_p1 = center_a .- p1
            # a_p2 = center_a .- p2
            # x = LA.dot(a_p1, p2_p1) * LA.dot(a_p2, p2_p1)
            # return x < zero(x)
        # end
    # end
# end

# is_colliding(a::GB.Line{N}, b::GB.HyperSphere{N}) where {N} = is_colliding(b, a)

#####
# HyperSphere vs. HyperSphere
#####

is_colliding(a::GB.HyperSphere, b::GB.HyperSphere, position_ba, direction_ba) = LA.dot(position_ba, position_ba) <= (a.r + b.r) ^ 2

# function is_colliding(a::GB.HyperSphere{N}, b::GB.HyperSphere{N}) where {N}
    # center_a = get_center(a)
    # center_b = get_center(b)
    # ba = center_b .- center_a
    # return LA.dot(ba, ba) <= (a.r + b.r) ^ 2
# end

#####
# HyperRectangle vs. Point
#####

is_colliding(a::GB.Rect, b::GB.Vec, position_ba, direction_ba) = minimum(a) <= position_ba <= maximum(a)

is_colliding(a::GB.Vec, b::GB.Rect, position_ba, direction_ba) = is_colliding(b, a, -position_ba, get_reverse_direction(direction_ba))

# is_colliding(a::GB.HyperRectangle{N}, b::GB.Vec{N}) where {N} = minimum(a) <= b <= maximum(a)

# is_colliding(a::GB.Vec{N}, b::GB.HyperRectangle{N}) where {N} = is_colliding(b, a)

#####
# HyperRectangle vs. Line
#####

function is_colliding(a::GB.Rect, b::GB.Line, position_ba, direction_ba) = minimum(a) <= position_ba <= maximum(a)
    if is_colliding(a, get_point(b, 1)) || is_colliding(a, get_point(b, 2))
        return true
    else
        half_widths = get_half_widths(a)
        half_width = half_widths[1]
        half_height = half_widths[2]

        PointType = typeof(b.points[1])
        VecType = typeof(position_ba)
        x = position_ba[1]
        y = position_ba[2]

        horizontal_line = GB.Line(PointType(-half_width, 0), PointType(half_width, 0))
        vertical_line = GB.Line(PointType(-half_height, 0), PointType(half_height, 0))

        direction_ba_minus_90_rotated = get_minus_90_rotated(direction_ba)
        return any([is_colliding(horizontal_line, b, VecType(x, y + half_height), direction_ba),
                    is_colliding(horizontal_line, b, VecType(x, y - half_height), direction_ba),
                    is_colliding(vertical_line, b, VecType(x + half_width, y), direction_ba_minus_90_rotated),
                    is_colliding(vertical_line, b, VecType(x - half_width, y), direction_ba_minus_90_rotated)])
    end
end

is_colliding(a::GB.Line, b::GB.Rect, position_ba, direction_ba) = is_colliding(b, a, -position_ba, get_reverse_direction(direction_ba))

# function is_colliding(a::GB.HyperRectangle{N}, b::GB.Line{N}) where {N}
    # if is_colliding(a, get_point(b, 1)) || is_colliding(a, get_point(b, 2))
        # return true
    # else
        # lines = get_lines(a)
        # return any(line -> is_colliding(line, b), lines)
    # end
# end

# is_colliding(a::GB.Line{N}, b::GB.HyperRectangle{N}) where {N} = is_colliding(b, a)

#####
# HyperRectangle vs. HyperSphere
#####

function is_colliding(a::GB.Rect, b::GB.HyperSphere, position_ba, direction_ba)
    if !(minimum(a) .- b.r <= position_ba <= maximum(a) .+ b.r)
        return false
    else
        bottom_left = get_bottom_left(a)
        position_bc = position_ba .- bottom_left
        if x < bottom_left[1] && y < bottom_left[2] && LA.dot(position_bc, position_bc) > r_squared
            return false
        end

        bottom_right = get_bottom_right(a)
        position_bc = position_ba .- bottom_right
        if x > bottom_left[1] && y < bottom_left[2] && LA.dot(position_bc, position_bc) > r_squared
            return false
        end

        top_right = get_top_right(a)
        position_bc = position_ba .- top_right
        if x > bottom_left[1] && y > bottom_left[2] && LA.dot(position_bc, position_bc) > r_squared
            return false
        end

        top_left = get_top_left(a)
        position_bc = position_ba .- top_left
        if x < bottom_left[1] && y > bottom_left[2] && LA.dot(position_bc, position_bc) > r_squared
            return false
        end

        return true
    end
end

is_colliding(a::GB.HyperSphere, b::GB.Rect, position_ba, direction_ba) = is_colliding(b, a, -position_ba, get_reverse_direction(direction_ba))

# function is_colliding(a::GB.HyperRectangle{N}, b::GB.HyperSphere{N}) where {N}
    # center_a = get_center(a)
    # center_b = get_center(b)
    # ba = center_b .- center_a

    # half_widths = get_half_widths(a)
    # ba_clamped = clamp.(ba, -half_widths, half_widths)

    # closest_point = center_a .+ ba_clamped

    # return is_colliding(b, closest_point)
# end

# is_colliding(a::GB.HyperSphere{N}, b::GB.HyperRectangle{N}) where {N} = is_colliding(b, a)

#####
# HyperRectangle vs. HyperRectangle
#####

function is_colliding(a::GB.Rect, b::GB.Rect, position_ba, direction_ba)
    half_widths_b = get_half_widths(b)
    normal_direction_ba = get_plus_90_rotated(direction_ba)
    half_width_b = half_widths[1]
    half_height_b = half_widths[2]

    bottom_left = position_ba .- half_width_b .* direction_ba .- half_height_b .* normal_direction_ba
    bottom_right = position_ba .+ half_width_b .* direction_ba .- half_height_b .* normal_direction_ba
    top_right = position_ba .+ half_width_b .* direction_ba .+ half_height_b .* normal_direction_ba
    top_left = position_ba .- half_width_b .* direction_ba .+ half_height_b .* normal_direction_ba

    half_widths_a = get_half_widths(a)
    half_width_a = half_widths_a[1]
    half_height_a = half_widths_a[2]

    min_x_b, max_x_b = extrema((bottom_left[1], bottom_right[1], top_right[1], top_left[1]))
    min_y_b, max_y_b = extrema((bottom_left[2], bottom_right[2], top_right[2], top_left[2]))

    return !((half_width_a < min_x_b) || (max_x_b < -half_width_a) || (half_height_a < min_y_b) || (max_y_b < -half_height_a))
end

# function is_colliding(a::GB.HyperRectangle{N}, b::GB.HyperRectangle{N}) where {N}
    # for i in 1:N
        # if (a.origin[i] + a.widths[i] < b.origin[i]) || (b.origin[i] + b.widths[i] < a.origin[i])
            # return false
        # end
    # end

    # return true
# end
