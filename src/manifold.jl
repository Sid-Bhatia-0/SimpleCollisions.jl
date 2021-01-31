const SQRT_2 = sqrt(2)

struct Manifold{T}
    penetration::T
    normal::GB.Vec2{T}
    contact::GB.Vec2{T}
end

get_penetration(manifold::Manifold) = manifold.penetration
get_normal(manifold::Manifold) = manifold.normal
get_contact(manifold::Manifold) = manifold.contact

function Manifold(a::RigidBody, b::RigidBody)
    shape_a = get_shape(a)
    shape_b = get_shape(b)

    position_a = get_position(a)
    position_b = get_position(b)
    position_ba = position_b .- position_a

    direction_a = get_direction(a)
    direction_b = get_direction(b)
    direction_ba = get_relative_direction(direction_a, direction_b)

    return Manifold(shape_a, shape_b, position_ba, direction_ba)
end

#####
# HyperSphere vs. HyperSphere
#####

function Manifold(a::GB.HyperSphere, b::GB.HyperSphere, position_ba, direction_ba)
    if all(ab -> isapprox(ab[1], ab[2], atol = ATOL), zip(zero(position_ba), position_ba))
        penetration = a.r + b.r
        normal = GB.unit(typeof(position_ba), 1)
        contact = zero(position_ba)
        return Manifold(penetration, normal, contact)
    else
        norm_ba = LA.norm(position_ba)
        penetration = a.r + b.r - norm_ba
        normal = position_ba ./ norm_ba
        contact = (a.r - norm_ba / 2) .* normal
        return Manifold(penetration, normal, contact)
    end
end

# function Manifold(a::GB.HyperSphere{N}, b::GB.HyperSphere{N}) where {N}
    # center_a = get_center(a)
    # center_b = get_center(b)
    # ba = center_b .- center_a

    # if all(ab -> isapprox(ab[1], ab[2], atol = ATOL), zip(center_a, center_b))
        # penetration = a.r + b.r
        # return Manifold(penetration, GB.unit(typeof(ba), 1), center_a)
    # else
        # norm_ba = LA.norm(ba)
        # penetration = a.r + b.r - norm_ba
        # normal = ba ./ norm_ba
        # return Manifold(penetration, normal, center_a .+ (a.r - norm_ba / 2) .* normal)
    # end
# end

#####
# HyperRectangle vs. HyperSphere
#####

function get_region(a, b, c)
    if c < a
        return 1
    elseif c > b
        return 3
    else
        return 2
    end
end

function get_penetration(a::GB.HyperRectangle, b::GB.Vec)
    bottom_left = minimum(a)
    top_right = maximum(a)
    min_penetration = Inf
    min_penetration_dim = -1
    min_penetration_dir = 0

    for i in 1:N
        penetration, dir = findmin((b[i] - bottom_left[i], top_right[i] - b[i]))
        if penetration < min_penetration
            min_penetration = penetration
            min_penetration_dim = i
            min_penetration_dir = dir == 1 ? -1 : 1
        end
    end

    return (min_penetration, min_penetration_dim, min_penetration_dir)
end

function Manifold(a::GB.Rect, b::GB.HyperSphere, position_ba, direction_ba)
    # center_b = get_center(b)
    bottom_left = minimum(a)
    top_right = maximum(a)
    region = get_region.(bottom_left, top_right, position_ba)

    if region === GB.Vec(1, 2)
        penetration = b.r - (bottom_left[1] - position_ba[1])
        normal = -GB.unit(typeof(position_ba), 1)
        contact = position_ba .+ (b.r - penetration / 2) .* -normal
        return Manifold(penetration, normal, contact)

    elseif region === GB.Vec(2, 1)
        penetration = b.r - (bottom_left[2] - position_ba[2])
        normal = -GB.unit(typeof(position_ba), 2)
        contact = position_ba .+ (b.r - penetration / 2) .* -normal
        return Manifold(penetration, normal, contact)

    elseif region === GB.Vec(2, 3)
        penetration = b.r - (position_ba[2] - top_right[2])
        normal = GB.unit(typeof(position_ba), 2)
        contact = position_ba .+ (b.r - penetration / 2) .* -normal
        return Manifold(penetration, normal, contact)

    elseif region === GB.Vec(3, 2)
        penetration = b.r - (position_ba[1] - top_right[1])
        normal = GB.unit(typeof(position_ba), 1)
        contact = position_ba .+ (b.r - penetration / 2) .* -normal
        return Manifold(penetration, normal, contact)

    elseif region === GB.Vec(1, 1)
        closest_point = bottom_left
        if all(ab -> isapprox(ab[1], ab[2], atol = ATOL), zip(position_ba, closest_point))
            penetration = b.r
            normal = typeof(position_ba)(-SQRT_2, -SQRT_2)
            contact = closest_point
            return Manifold(penetration, normal, contact)
        else
            vec = position_ba .- closest_point
            d = LA.norm(vec) # this won't be zero coz we already checked for it
            penetration = b.r - d
            normal = vec ./ d
            contact = position_ba .+ (b.r - penetration / 2) .* -normal
            return Manifold(penetration, normal, contact)
        end

    elseif region === GB.Vec(1, 3)
        closest_point = get_top_left(a)
        if all(ab -> isapprox(ab[1], ab[2], atol = ATOL), zip(position_ba, closest_point))
            penetration = b.r
            normal = typeof(position_ba)(-SQRT_2, SQRT_2)
            contact = closest_point
            return Manifold(penetration, normal, contact)
        else
            vec = position_ba .- closest_point
            d = LA.norm(vec) # this won't be zero coz we already checked for it
            penetration = b.r - d
            normal = vec ./ d
            contact = position_ba .+ (b.r - penetration / 2) .* -normal
            return Manifold(penetration, normal, contact)
        end

    elseif region === GB.Vec(3, 1)
        closest_point = get_bottom_right(a)
        if all(ab -> isapprox(ab[1], ab[2], atol = ATOL), zip(position_ba, closest_point))
            penetration = b.r
            normal = typeof(position_ba)(SQRT_2, -SQRT_2)
            contact = closest_point
            return Manifold(penetration, normal, contact)
        else
            vec = position_ba .- closest_point
            d = LA.norm(vec) # this won't be zero coz we already checked for it
            penetration = b.r - d
            normal = vec ./ d
            contact = position_ba .+ (b.r - penetration / 2) .* -normal
            return Manifold(penetration, normal, contact)
        end

    elseif region === GB.Vec(3, 3)
        closest_point = top_right
        if all(ab -> isapprox(ab[1], ab[2], atol = ATOL), zip(position_ba, closest_point))
            penetration = b.r
            normal = typeof(position_ba)(SQRT_2, SQRT_2)
            contact = closest_point
            return Manifold(penetration, normal, contact)
        else
            vec = position_ba .- closest_point
            d = LA.norm(vec) # this won't be zero coz we already checked for it
            penetration = b.r - d
            normal = vec ./ d
            contact = position_ba .+ (b.r - penetration / 2) .* -normal
            return Manifold(penetration, normal, contact)
        end

    elseif region === GB.Vec(2, 2)
        penetration, dim, dir = get_penetration(a, position_ba)
        penetration = penetration + b.r
        normal = dir * GB.unit(typeof(position_ba), dim)
        contact = position_ba .+ (b.r - penetration / 2) * -normal
        return Manifold(penetration, normal, contact)
    end
end

function Manifold(a::GB.HyperSphere, b::GB.Rect, position_ba, direction_ba)
    manifold = Manifold(b, a, -position_ba, get_reverse_direction(direction_ba))
    Manifold(manifold.penetration, -manifold.normal, manifold.contact)
end

# function Manifold(a::GB.HyperRectangle{2}, b::GB.HyperSphere{2})
    # center_b = get_center(b)
    # bottom_left = minimum(a)
    # top_right = maximum(a)
    # region = get_region.(bottom_left, top_right, center_b)

    # if region === GB.Vec(1, 2)
        # penetration = b.r - (bottom_left[1] - center_b[1])
        # normal = -GB.unit(typeof(center_b), 1)
        # contact = center_b .+ (b.r - penetration / 2) .* -normal
        # return Manifold(penetration, normal, contact)

    # elseif region === GB.Vec(2, 1)
        # penetration = b.r - (bottom_left[2] - center_b[2])
        # normal = -GB.unit(typeof(center_b), 2)
        # contact = center_b .+ (b.r - penetration / 2) .* -normal
        # return Manifold(penetration, normal, contact)

    # elseif region === GB.Vec(2, 3)
        # penetration = b.r - (center_b[2] - top_right[2])
        # normal = GB.unit(typeof(center_b), 2)
        # contact = center_b .+ (b.r - penetration / 2) .* -normal
        # return Manifold(penetration, normal, contact)

    # elseif region === GB.Vec(3, 2)
        # penetration = b.r - (center_b[1] - top_right[1])
        # normal = GB.unit(typeof(center_b), 1)
        # contact = center_b .+ (b.r - penetration / 2) .* -normal
        # return Manifold(penetration, normal, contact)

    # elseif region === GB.Vec(1, 1)
        # closest_point = bottom_left
        # if all(ab -> isapprox(ab[1], ab[2], atol = ATOL), zip(center_b, closest_point))
            # penetration = b.r
            # normal = typeof(center_b)(-SQRT_2, -SQRT_2)
            # contact = closest_point
            # return Manifold(penetration, normal, contact)
        # else
            # vec = center_b .- closest_point
            # d = LA.norm(vec) # this won't be zero coz we already checked for it
            # penetration = b.r - d
            # normal = vec ./ d
            # contact = center_b .+ (b.r - penetration / 2) .* -normal
            # return Manifold(penetration, normal, contact)
        # end

    # elseif region === GB.Vec(1, 3)
        # closest_point = get_top_left(a)
        # if all(ab -> isapprox(ab[1], ab[2], atol = ATOL), zip(center_b, closest_point))
            # penetration = b.r
            # normal = typeof(center_b)(-SQRT_2, SQRT_2)
            # contact = closest_point
            # return Manifold(penetration, normal, contact)
        # else
            # vec = center_b .- closest_point
            # d = LA.norm(vec) # this won't be zero coz we already checked for it
            # penetration = b.r - d
            # normal = vec ./ d
            # contact = center_b .+ (b.r - penetration / 2) .* -normal
            # return Manifold(penetration, normal, contact)
        # end

    # elseif region === GB.Vec(3, 1)
        # closest_point = get_bottom_right(a)
        # if all(ab -> isapprox(ab[1], ab[2], atol = ATOL), zip(center_b, closest_point))
            # penetration = b.r
            # normal = typeof(center_b)(SQRT_2, -SQRT_2)
            # contact = closest_point
            # return Manifold(penetration, normal, contact)
        # else
            # vec = center_b .- closest_point
            # d = LA.norm(vec) # this won't be zero coz we already checked for it
            # penetration = b.r - d
            # normal = vec ./ d
            # contact = center_b .+ (b.r - penetration / 2) .* -normal
            # return Manifold(penetration, normal, contact)
        # end

    # elseif region === GB.Vec(3, 3)
        # closest_point = get_bottom_right(a)
        # if all(ab -> isapprox(ab[1], ab[2], atol = ATOL), zip(center_b, closest_point))
            # penetration = b.r
            # normal = typeof(center_b)(SQRT_2, SQRT_2)
            # contact = closest_point
            # return Manifold(penetration, normal, contact)
        # else
            # vec = center_b .- closest_point
            # d = LA.norm(vec) # this won't be zero coz we already checked for it
            # penetration = b.r - d
            # normal = vec ./ d
            # contact = center_b .+ (b.r - penetration / 2) .* -normal
            # return Manifold(penetration, normal, contact)
        # end

    # elseif region === GB.Vec(2, 2)
        # penetration, dim, dir = get_penetration(a, center_b)
        # penetration = penetration + b.r
        # normal = dir * GB.unit(typeof(center_b), dim)
        # contact = center_b .+ (b.r - penetration / 2) * -normal
        # return Manifold(penetration, normal, contact)
    # end
# end

# function Manifold(a::GB.HyperSphere{N}, b::GB.HyperRectangle{N}) where {N}
    # manifold = Manifold(b, a)
    # Manifold(manifold.penetration, -manifold.normal, manifold.contact)
# end

#####
# HyperRectangle vs. HyperRectangle
#####

function Manifold(a::GB.HyperRectangle{N}, b::GB.HyperRectangle{N}) where {N}
    intersection = GB.intersect(a, b)
    penetration, dim = findmin(intersection.widths)
    contact = get_center(intersection)

    normal = GB.unit(typeof(a.origin), dim)
    if a.origin[dim] > b.origin[dim]
        normal = -normal
    end
    
    return Manifold(penetration, normal, contact)
end
