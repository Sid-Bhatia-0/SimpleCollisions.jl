#####
# Point vs. Point
#####

is_colliding(a::GB.Point{N}, b::GB.Point{N}) where {N} = a == b

#####
# Line vs. Point
#####

function is_colliding(a::GB.Line{N}, b::GB.Point{N}) where {N}
    p1 = a.points[1]
    p2 = a.points[2]
    p1_p2 = p1 .- p2
    b_p1 = p1 .- b
    b_p2 = p2 .- b
    d1_squared = LA.dot(b_p1, b_p1)
    d2_squared = LA.dot(b_p2, b_p2)
    d3_squared = LA.dot(p1_p2, p1_p2)
    return (d1_squared + d2_squared - d3_squared) ^ 2 ≈ 4 * d1_squared * d2_squared
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
