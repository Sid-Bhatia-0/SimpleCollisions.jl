#####
# Point vs. Point
#####

is_colliding(a::GB.Point{N}, b::GB.Point{N}) where {N} = a == b

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
# HyperRectangle vs. HyperSphere
#####

function is_colliding(a::GB.HyperRectangle{N}, b::GB.HyperSphere{N}) where {N}
    a_center = get_center(a)
    ba = b.center .- a_center

    half_widths = get_half_widths(a)
    ba_clamped = clamp(ba, -half_widths, half_widths)

    closest_point = GB.Point(a_center .+ ba_clamped)

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
