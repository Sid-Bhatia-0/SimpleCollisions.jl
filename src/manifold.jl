struct Manifold{T}
    penetration::T
    normal::GB.Vec2{T}
end

#####
# HyperSphere vs. HyperSphere
#####

function Manifold(a::GB.HyperSphere{N}, b::GB.HyperSphere{N}) where {N}
    ba = GB.Vec(b.center .- a.center)
    d = LA.norm(ba)
    penetration = a.r + b.r - d

    if d == zero(d)
        return Manifold(penetration, GB.unit(typeof(ba), 1))
    else
        return Manifold(penetration, LA.normalize(ba))
    end
end

#####
# HyperRectangle vs. HyperSphere
#####

function Manifold(a::GB.HyperRectangle{N}, b::GB.HyperSphere{N}) where {N}
    center_a = get_center(a)
    center_b = get_center(b)
    ba = GB.Vec(center_b .- center_a)

    half_widths = get_half_widths(a)
    ba_clamped = clamp.(ba, -half_widths, half_widths)

    d = LA.norm(ba_clamped)
    closest_point = GB.Point(center_a .+ ba_clamped)
    normal = LA.norm(GB.Vec(center_b .- closest_point))
    penetration = b.r - LA.norm(GB.Vec(closest_point .- center_b))

    if d == zero(d)
        return Manifold(penetration, GB.unit(typeof(ba), 1))
    else
        return Manifold(penetration, LA.normalize(ba_clamped))
    end
end

function Manifold(a::GB.HyperSphere{N}, b::GB.HyperRectangle{N}) where {N}
    manifold = Manifold(b, a)
    Manifold(manifold.penetration, -manifold.normal)
end

#####
# HyperRectangle vs. HyperRectangle
#####

function Manifold(a::GB.HyperRectangle{N}, b::GB.HyperRectangle{N}) where {N}
    intersection = GB.intersect(a, b)
    penetration, dim = findmin(intersection.widths)

    if a.origin[dim] <= b.origin[dim]
        return Manifold(penetration, GB.unit(typeof(a.origin), dim))
    else
        return Manifold(penetration, -GB.unit(typeof(a.origin), dim))
    end
end
