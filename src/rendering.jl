get_plotting_shape(circle::StdCircle{T}, pos::GB.Vec{2, T}) where {T} = GB.Circle(convert(GB.Point, pos), get_radius(circle))

get_plotting_shape(rect::StdRect{T}, pos::GB.Vec{2, T}) where {T} = GB.Rect(pos .+ get_bottom_left(rect)..., get_width(rect), get_height(rect))
