#####
# Point vs. Point
#####

is_colliding(a::GB.Point{N}, b::GB.Point{N}) where {N} = a == b

#####
# Line vs. Point
#####

function is_colliding(a::GB.Line{N}, b::GB.Point{N}) where {N}
    T = eltype(b)
    p1 = a.points[1]
    p2 = a.points[2]
    b_p1 = p1 .- b
    b_p2 = p2 .- b
    return LA.dot(b_p1, b_p2) <= zero(T) && get_area(b, p1, p2) ≈ zero(T)
end

is_colliding(a::GB.Point{N}, b::GB.Line{N}) where {N} = is_colliding(b, a)

#####
# HyperSphere vs. Point
#####

function is_colliding(a::GB.HyperSphere{N}, b::GB.Point{N}) where {N}
    ba = b .- a.center
    return LA.dot(ba, ba) <= a.r ^ 2
end

is_colliding(a::GB.Point{N}, b::GB.HyperSphere{N}) where {N} = is_colliding(b, a)

#####
# HyperSphere vs. Line
#####

function is_colliding(a::GB.HyperSphere{N}, b::GB.Line{N}) where {N}
    p1 = b.points[1]
    p2 = b.points[2]
    if is_colliding(a, p1) || is_colliding(a, p2)
        return true
    else
        center_a = get_center(a)
        area = get_area(center_a, p1, p2)
        p1_p2 = p2 .- p1
        height_squared = 4 * area ^ 2 / LA.dot(p1_p2, p1_p2)
        if height_squared > a.r ^ 2
            return false
        else
            p1_a = center_a .- p1
            p2_a = center_a .- p2
            return LA.dot(p1_a, p1_p2) * LA.dot(p2_a, p1_p2) < zero(eltype(p1_a))
        end
    end
end

is_colliding(a::GB.Line{N}, b::GB.HyperSphere{N}) where {N} = is_colliding(b, a)

#####
# HyperSphere vs. HyperSphere
#####

function is_colliding(a::GB.HyperSphere{N}, b::GB.HyperSphere{N}) where {N}
    ba = b.center .- a.center
    return LA.dot(ba, ba) <= (a.r + b.r) ^ 2
end

#####
# HyperRectangle vs. Point
#####

is_colliding(a::GB.HyperRectangle{N}, b::GB.Point{N}) where {N} = minimum(a) <= b <= maximum(a)

is_colliding(a::GB.Point{N}, b::GB.HyperRectangle{N}) where {N} = is_colliding(b, a)

#####
# HyperRectangle vs. HyperSphere
#####

function is_colliding(a::GB.HyperRectangle{N}, b::GB.HyperSphere{N}) where {N}
    center_a = get_center(a)
    center_b = get_center(b)
    ba = GB.Vec(center_b .- center_a)

    half_widths = get_half_widths(a)
    ba_clamped = clamp.(ba, -half_widths, half_widths)

    closest_point = GB.Point(center_a .+ ba_clamped)

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
