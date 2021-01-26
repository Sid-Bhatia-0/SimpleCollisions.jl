struct Manifold{T}
    penetration::T
    normal::GB.Vec2{T}
    contact::GB.Vec2{T}
end

get_penetration(manifold::Manifold) = manifold.penetration
get_normal(manifold::Manifold) = manifold.normal
get_contact(manifold::Manifold) = manifold.contact

#####
# HyperSphere vs. HyperSphere
#####

function Manifold(a::GB.HyperSphere{N}, b::GB.HyperSphere{N}) where {N}
    ba = GB.Vec(b.center .- a.center)
    d = LA.norm(ba)
    penetration = a.r + b.r - d

    # assuming ba is not the zero vector
    normal = ba ./ d

    return Manifold(penetration, normal, a.center .+ (a.r - d/2) * normal)
end

#####
# HyperRectangle vs. HyperSphere
#####

function Manifold(a::GB.HyperRectangle{N}, b::GB.HyperSphere{N}) where {N}
    center_a = get_center(a)
    center_b = get_center(b)
    # assuming ba is not the zero vector
    ba = GB.Vec(center_b .- center_a)

    if !is_colliding(a, center_b)
        half_widths = get_half_widths(a)
        ba_clamped = clamp.(ba, -half_widths, half_widths)

        closest_point = center_a .+ ba_clamped
        normal = GB.Vec(center_b .- closest_point)
        d = LA.norm(normal)
        penetration = b.r - d # max penetration possible is b.r since d >= 0
        normal = normal ./ d # assuming d is not zero
        contact = closest_point .+ (d / 2) * normal

        return Manifold(penetration, normal, contact)
    else
        penetration = b.r
        d = LA.norm(ba)
        normal = normal ./ d # assuming d is not zero
        contact = center_a .+ (d - b.r/2) * normal
        return Manifold(b.r, LA.normalize(ba), contact)
    end
end

function Manifold(a::GB.HyperSphere{N}, b::GB.HyperRectangle{N}) where {N}
    manifold = Manifold(b, a)
    Manifold(manifold.penetration, -manifold.normal, manifold.contact)
end

#####
# HyperRectangle vs. HyperRectangle
#####

function Manifold(a::GB.HyperRectangle{N}, b::GB.HyperRectangle{N}) where {N}
    intersection = GB.intersect(a, b)
    penetration, dim = findmin(intersection.widths)
    contact = get_center(intersection)

    normal = GB.unit(typeof(a.origin), dim)
    if a.origin[dim] <= b.origin[dim]
        normal = -normal
    end
    
    return Manifold(penetration, normal, contact)
end
