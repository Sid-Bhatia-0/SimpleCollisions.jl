#####
# StandardLine vs. StandardLine
#####

function is_colliding(a::StandardLine{T}, b::StandardLine{T}, pos_ba::Vector2D{T}) where {T}
    y = pos_ba[2]
    if iszero(y)
        x = pos_ba[1]
        half_length_a = get_half_length(a)
        half_length_b = get_half_length(b)
        return !((x + half_length_b <= -half_length_a) || (x - half_length_b >= half_length_a))
    else
        return false
    end
end

function is_colliding(l1::StandardLine{T}, l2::StandardLine{T}, pos::Vector2D{T}, dir::Vector2D{T}) where {T}
    half_length = get_half_length(l1)

    tail_l2, head_l2 = get_vertices(l2, pos, dir)

    x1 = tail_l2[1]
    y1 = tail_l2[2]

    x2 = head_l2[1]
    y2 = head_l2[2]

    t = y1 / (y1 - y2)

    if isfinite(t)
        x_intersection = x1 + t * (x2 - x1)
        return (zero(t) < t < one(t)) && (-half_length < x_intersection < half_length)
    else
        return (iszero(y1)) && (iszero(y2)) && !((half_length <= x1) || (x2 <= -half_length))
    end
end

#####
# StandardCircle vs. StandardPoint
#####

function is_inside(circle::StandardCircle{T}, pos::Vector2D{T}) where {T}
    radius = get_radius(circle)
    return LA.dot(pos, pos) < radius * radius
end

is_colliding(circle::StandardCircle{T}, point::StandardPoint{T}, pos::Vector2D{T}) where {T} = is_inside(circle, pos)
is_colliding(point::StandardPoint{T}, circle::StandardCircle{T}, pos::Vector2D{T}) where {T} = is_inside(circle, pos) # no need to reverse pos because of symmetry

is_colliding(circle::StandardCircle{T}, point::StandardPoint{T}, pos::Vector2D{T}, dir::Vector2D{T}) where {T} = is_inside(circle, pos)
is_colliding(point::StandardPoint{T}, circle::StandardCircle{T}, pos::Vector2D{T}, dir::Vector2D{T}) where {T} = is_inside(circle, pos) # no need to reverse pos because of symmetry

#####
# StandardCircle vs. StandardLine
#####

function get_projection(line::StandardLine{T}, pos::Vector2D{T}) where {T}
    half_length = get_half_length(line)
    x = pos[1]
    if x < -half_length
        return get_tail(line)
    elseif x > half_length
        return get_head(line)
    else
        return typeof(pos)(x, zero(x))
    end
end

function is_colliding(line::StandardLine{T}, circle::StandardCircle{T}, pos::Vector2D{T}) where {T}
    projection = get_projection(line, pos)
    vec = pos - projection
    radius = get_radius(circle)
    return LA.dot(vec, vec) < radius * radius
end

is_colliding(circle::StandardCircle{T}, line::StandardLine{T}, pos::Vector2D{T}) where {T} = is_colliding(line, circle, pos) # no need to reverse pos because of symmetry

is_colliding(line::StandardLine{T}, circle::StandardCircle{T}, pos::Vector2D{T}, dir::Vector2D{T}) where {T} = is_colliding(line, circle, pos)
is_colliding(circle::StandardCircle{T}, line::StandardLine{T}, pos::Vector2D{T}, dir::Vector2D{T}) where {T} = is_colliding(line, circle, invert(pos, dir)...)

#####
# StandardCircle vs. StandardCircle
#####

function is_colliding(c1::StandardCircle{T}, c2::StandardCircle{T}, pos::Vector2D{T}) where {T}
    r1 = get_radius(c1)
    r2 = get_radius(c2)
    r = r1 + r2
    return LA.dot(pos, pos) < r * r
end

is_colliding(c1::StandardCircle{T}, c2::StandardCircle{T}, pos::Vector2D{T}, dir::Vector2D{T}) where {T} = is_colliding(c1, c2, pos)

#####
# StandardRect vs. StandardPoint
#####

function is_inside(rect::StandardRect{T}, pos::Vector2D{T}) where {T}
    half_width = get_half_width(rect)
    half_height = get_half_height(rect)

    x = pos[1]
    y = pos[2]

    return (-half_width < x < half_width) && (-half_height < y < half_height)
end

is_colliding(rect::StandardRect{T}, point::StandardPoint{T}, pos::Vector2D{T}) where {T} = is_inside(rect, pos)
is_colliding(point::StandardPoint{T}, rect::StandardRect{T}, pos::Vector2D{T}) where {T} = is_inside(rect, pos) # no need to reverse pos because of symmetry

is_colliding(rect::StandardRect{T}, point::StandardPoint{T}, pos::Vector2D{T}, dir::Vector2D{T}) where {T} = is_inside(rect, pos)
is_colliding(point::StandardPoint{T}, rect::StandardRect{T}, pos::Vector2D{T}, dir::Vector2D{T}) where {T} = is_colliding(rect, point, invert(pos, dir)...)

#####
# StandardRect vs. StandardLine
#####

function separating_axis_exists(rect::StandardRect{T}, line::StandardLine{T}, pos::Vector2D{T}) where {T}
    half_height = get_half_height(rect)
    half_width = get_half_width(rect)

    x = pos[1]
    y = pos[2]
    half_length = get_half_length(line)

    return (y <= -half_height) || (y >= half_height) || (x + half_length <= -half_width) || (x - half_length >= half_width)
end

is_colliding(rect::StandardRect{T}, line::StandardLine{T}, pos::Vector2D{T}) where {T} = !separating_axis_exists(rect, line, pos)
is_colliding(line::StandardLine{T}, rect::StandardRect{T}, pos::Vector2D{T}) where {T} = is_colliding(rect, line, pos) # no need to reverse pos because of symmetry

function separating_axis_exists(rect::StandardRect{T}, line::StandardLine{T}, pos::Vector2D{T}, dir::Vector2D{T}) where {T}
    half_width = get_half_width(rect)
    half_height = get_half_height(rect)

    tail, head = get_vertices(line, pos, dir)

    min_x, max_x = minmax(tail[1], head[1])
    min_y, max_y = minmax(tail[2], head[2])

    return (max_x <= -half_width) || (min_x >= half_width) || (max_y <= -half_height) || (min_y >= half_height)
end

function separating_axis_exists(line::StandardLine{T}, rect::StandardRect{T}, pos::Vector2D{T}, dir::Vector2D{T}) where {T}
    half_length = get_half_length(line)

    bottom_left, bottom_right, top_right, top_left = get_vertices(rect, pos, dir)

    min_x, max_x = extrema((bottom_left[1], bottom_right[1], top_right[1], top_left[1]))
    min_y, max_y = extrema((bottom_left[2], bottom_right[2], top_right[2], top_left[2]))

    return (max_x <= -half_length) || (min_x >= half_length) || (max_y <= zero(T)) || (min_y >= zero(T))
end

is_colliding(rect::StandardRect{T}, line::StandardLine{T}, pos::Vector2D{T}, dir::Vector2D{T}) where {T} = !(separating_axis_exists(rect, line, pos, dir) || separating_axis_exists(line, rect, invert(pos, dir)...))
is_colliding(line::StandardLine{T}, rect::StandardRect{T}, pos::Vector2D{T}, dir::Vector2D{T}) where {T} = !(separating_axis_exists(line, rect, pos, dir) || separating_axis_exists(rect, line, invert(pos, dir)...))

#####
# StandardRect vs. StandardCircle
#####

get_projection(rect::StandardRect{T}, pos::Vector2D{T}) where {T} = clamp.(pos, get_bottom_left(rect), get_top_right(rect))

function is_colliding(rect::StandardRect{T}, circle::StandardCircle{T}, pos::Vector2D{T}) where {T}
    projection = get_projection(rect, pos)
    vec = pos - projection
    radius = get_radius(circle)
    return LA.dot(vec, vec) < radius * radius
end

is_colliding(circle::StandardCircle{T}, rect::StandardRect{T}, pos::Vector2D{T}) where {T} = is_colliding(rect, circle, pos) # no need to invert pos because of symmetry

is_colliding(rect::StandardRect{T}, circle::StandardCircle{T}, pos::Vector2D{T}, dir::Vector2D{T}) where {T} = is_colliding(rect, circle, pos)
is_colliding(circle::StandardCircle{T}, rect::StandardRect{T}, pos::Vector2D{T}, dir::Vector2D{T}) where {T} = is_colliding(rect, circle, invert(pos, dir)...)

#####
# StandardRect vs. StandardRect
#####

function separating_axis_exists(r1::StandardRect{T}, r2::StandardRect{T}, pos::Vector2D{T}) where {T}
    half_width_r1 = get_half_width(r1)
    half_height_r1 = get_half_height(r1)

    half_width_r2 = get_half_width(r2)
    half_height_r2 = get_half_height(r2)

    x = pos[1]
    y = pos[2]

    return (x + half_width_r2 <= -half_width_r1) || (x - half_width_r2 >= half_width_r1) || (y + half_height_r2 <= -half_height_r1) || (y - half_height_r2 >= half_height_r1)
end

is_colliding(a::StandardRect{T}, b::StandardRect{T}, pos_ba::Vector2D{T}) where {T} = !separating_axis_exists(a, b, pos_ba)

function separating_axis_exists(r1::StandardRect{T}, r2::StandardRect{T}, pos::Vector2D{T}, dir::Vector2D{T}) where {T}
    half_width_r1 = get_half_width(r1)
    half_height_r1 = get_half_height(r1)

    bottom_left_r2, bottom_right_r2, top_right_r2, top_left_r2 = get_vertices(r2, pos, dir)

    min_x_r2, max_x_r2 = extrema((bottom_left_r2[1], bottom_right_r2[1], top_right_r2[1], top_left_r2[1]))
    min_y_r2, max_y_r2 = extrema((bottom_left_r2[2], bottom_right_r2[2], top_right_r2[2], top_left_r2[2]))

    return ((half_width_r1 <= min_x_r2) || (max_x_r2 <= -half_width_r1) || (half_height_r1 <= min_y_r2) || (max_y_r2 <= -half_height_r1))
end

is_colliding(r1::StandardRect{T}, r2::StandardRect{T}, pos::Vector2D{T}, dir::Vector2D{T}) where {T} = !(separating_axis_exists(r1, r2, pos, dir) || separating_axis_exists(r2, r1, invert(pos, dir)...))
