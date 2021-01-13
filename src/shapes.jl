GB.area(a::GB.Rect2D) = prod(a.widths)
GB.area(a::GB.Circle) = π * a.r * a.r

get_mass(density::T, shape::GB.GeometryPrimitive{2}) where {T<:Number} = convert(T, density * GB.area(shape))

get_inertia(density::T, shape::GB.Circle) where {T<:Number} = convert(T, 0.5 * get_mass(density, shape) * shape.r * shape.r)
get_inertia(density::T, shape::GB.Rect2D) where {T<:Number} = convert(T, get_mass(density, shape) * LA.dot(shape.widths, shape.widths) / 12)
