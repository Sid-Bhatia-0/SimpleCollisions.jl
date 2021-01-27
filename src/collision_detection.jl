#####
# Point vs. Point
#####

is_colliding(a::GB.Vec{N}, b::GB.Vec{N}) where {N} = a ≈ b

#####
# Line vs. Point
#####

function is_colliding(a::GB.Line{N}, b::GB.Vec{N}) where {N}
    T = eltype(b)
    p1 = a.points[1]
    p2 = a.points[2]
    p1_b = p1 .- b
    p2_b = p2 .- b
    return LA.dot(p1_b, p2_b) <= zero(T) && get_area(b, p1, p2) ≈ zero(T)
end

is_colliding(a::GB.Vec{N}, b::GB.Line{N}) where {N} = is_colliding(b, a)

#####
# Line vs. Line
#####

get_point(a::GB.Line, i) = convert(GB.Vec, a.points[i])

function is_colliding(a::GB.Line{2}, b::GB.Line{2})
    p1 = get_point(a, 1)
    p2 = get_point(a, 2)
    q1 = get_point(b, 1)
    q2 = get_point(b, 2)
    p2_p1 = p2 .- p1
    q2_q1 = q2 .- q1

    normal_a = GB.Vec(-p2_p1[2], p2_p1[1])
    x = LA.dot(q2_q1, normal_a)
    if x ≈ zero(x)
        return is_colliding(a, q1) || is_colliding(a, q2)
    else
        normal_b = GB.Vec(-q2_q1[2], q2_q1[1])
        y = LA.dot(p2_p1, normal_b)
        g = LA.dot(q1 .- p1, normal_b) / y # note that since the lines are not parallel (since x is non-zero), y will not be zero
        h = LA.dot(p1 .- q1, normal_a) / x
        if zero(g) <= g <= one(g) && zero(h) <= h <= one(h)
            return true
        else
            return false
        end
    end
end

#####
# HyperSphere vs. Point
#####

function is_colliding(a::GB.HyperSphere{N}, b::GB.Vec{N}) where {N}
    center_a = get_center(a)
    ba = b .- center_a
    return LA.dot(ba, ba) <= a.r ^ 2
end

is_colliding(a::GB.Vec{N}, b::GB.HyperSphere{N}) where {N} = is_colliding(b, a)

#####
# HyperSphere vs. Line
#####

function is_colliding(a::GB.HyperSphere{N}, b::GB.Line{N}) where {N}
    p1 = get_point(b, 1)
    p2 = get_point(b, 2)
    if is_colliding(a, p1) || is_colliding(a, p2)
        return true
    else
        center_a = get_center(a)
        area = get_area(center_a, p1, p2)
        p2_p1 = p2 .- p1
        height_squared = 4 * area ^ 2 / LA.dot(p2_p1, p2_p1)
        if height_squared > a.r ^ 2
            return false
        else
            a_p1 = center_a .- p1
            a_p2 = center_a .- p2
            x = LA.dot(a_p1, p2_p1) * LA.dot(a_p2, p2_p1)
            return x < zero(x)
        end
    end
end

is_colliding(a::GB.Line{N}, b::GB.HyperSphere{N}) where {N} = is_colliding(b, a)

#####
# HyperSphere vs. HyperSphere
#####

function is_colliding(a::GB.HyperSphere{N}, b::GB.HyperSphere{N}) where {N}
    center_a = get_center(a)
    center_b = get_center(b)
    ba = center_b .- center_a
    return LA.dot(ba, ba) <= (a.r + b.r) ^ 2
end

#####
# HyperRectangle vs. Point
#####

is_colliding(a::GB.HyperRectangle{N}, b::GB.Vec{N}) where {N} = minimum(a) <= b <= maximum(a)

is_colliding(a::GB.Vec{N}, b::GB.HyperRectangle{N}) where {N} = is_colliding(b, a)

#####
# HyperRectangle vs. Line
#####

function is_colliding(a::GB.HyperRectangle{N}, b::GB.Line{N}) where {N}
    if is_colliding(a, get_point(b, 1)) || is_colliding(a, get_point(b, 2))
        return true
    else
        lines = get_lines(a)
        return any(line -> is_colliding(line, b), lines)
    end
end

is_colliding(a::GB.Line{N}, b::GB.HyperRectangle{N}) where {N} = is_colliding(b, a)

#####
# HyperRectangle vs. HyperSphere
#####

function is_colliding(a::GB.HyperRectangle{N}, b::GB.HyperSphere{N}) where {N}
    center_a = get_center(a)
    center_b = get_center(b)
    ba = center_b .- center_a

    half_widths = get_half_widths(a)
    ba_clamped = clamp.(ba, -half_widths, half_widths)

    closest_point = center_a .+ ba_clamped

    return is_colliding(b, closest_point)
end

is_colliding(a::GB.HyperSphere{N}, b::GB.HyperRectangle{N}) where {N} = is_colliding(b, a)

#####
# HyperRectangle vs. HyperRectangle
#####

function is_colliding(a::GB.HyperRectangle{N}, b::GB.HyperRectangle{N}) where {N}
    for i in 1:N
        if (a.origin[i] + a.widths[i] < b.origin[i]) || (b.origin[i] + b.widths[i] < a.origin[i])
            return false
        end
    end

    return true
end
