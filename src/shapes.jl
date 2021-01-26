get_center(a::GB.HyperSphere) = a.center
get_center(a::GB.HyperRectangle) = GB.Point(a.origin .+ a.widths ./ 2)
get_half_widths(a::GB.HyperRectangle) = a.widths ./ 2

GB.area(a::GB.Rect2D) = prod(a.widths)
GB.area(a::GB.Circle) = π * a.r * a.r
function get_area(p1::GB.Point2, p2::GB.Point2, p3::GB.Point2)
    x1 = p1[1]
    y1 = p1[2]
    x2 = p2[1]
    y2 = p2[2]
    x3 = p3[1]
    y3 = p3[2]
    return abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2
end

function get_lines(a::GB.Rect{2, T}) where {T}
    bottom_left = GB.Point(minimum(a))
    bottom_right = GB.Point(bottom_left .+ GB.Vec(a.widths[1], zero(T)))
    top_right = GB.Point(maximum(a))
    top_left = GB.Point(bottom_left .+ GB.Vec(zero(T), a.widths[2]))

    l1 = GB.Line(bottom_left, bottom_right)
    l2 = GB.Line(bottom_right, top_right)
    l3 = GB.Line(top_right, top_left)
    l4 = GB.Line(top_left, bottom_left)
    return (l1, l2, l3, l4)
end

get_mass(density::Number, shape::GB.GeometryPrimitive{2}) = density * GB.area(shape)

get_inertia(density::Number, shape::GB.Circle) = get_mass(density, shape) * shape.r * shape.r / 2
get_inertia(density::Number, shape::GB.Rect2D) = get_mass(density, shape) * LA.dot(shape.widths, shape.widths) / 12
