get_center(a::GB.HyperSphere) = a.center
get_center(a::GB.HyperRectangle) = a.origin .+ a.widths ./ 2
get_half_widths(a::GB.HyperRectangle) = a.widths ./ 2

GB.area(a::GB.Rect2D) = prod(a.widths)
GB.area(a::GB.Circle) = π * a.r * a.r

get_mass(density::Number, shape::GB.GeometryPrimitive{2}) = density * GB.area(shape)

get_inertia(density::Number, shape::GB.Circle) = get_mass(density, shape) * shape.r * shape.r / 2
get_inertia(density::Number, shape::GB.Rect2D) = get_mass(density, shape) * LA.dot(shape.widths, shape.widths) / 12
