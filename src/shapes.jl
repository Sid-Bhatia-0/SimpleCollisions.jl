get_point(a::GB.Line, i) = convert(GB.Vec, a.points[i])
get_center(a::GB.HyperSphere) = convert(GB.Vec, a.center)
get_center(a::GB.HyperRectangle) = convert(GB.Vec, a.origin .+ a.widths ./ 2)
get_half_widths(a::GB.HyperRectangle) = convert(GB.Vec, a.widths ./ 2)
get_bottom_left(a::GB.HyperRectangle) = minimum(a)
get_bottom_right(a::GB.HyperRectangle{N, T}) where {N, T} = minimum(a) .+ GB.Vec(GB.widths[1], zero(T))
get_top_right(a::GB.HyperRectangle) = maximum(a)
get_top_left(a::GB.HyperRectangle) = minimum(a) .+ GB.Vec(zero(T), GB.widths[2])

get_area(a::GB.Rect2D) = prod(a.widths)
get_area(a::GB.Circle) = π * a.r * a.r
function get_area(p1, p2, p3)
    x1 = p1[1]
    y1 = p1[2]
    x2 = p2[1]
    y2 = p2[2]
    x3 = p3[1]
    y3 = p3[2]
    return abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2
end

function get_edges(a::GB.Rect2D)
    bottom_left = convert(GB.Point, get_bottom_left(a))
    bottom_right = convert(GB.Point, get_bottom_right(a))
    top_right = convert(GB.Point, get_top_right(a))
    top_left = convert(GB.Point, get_top_left(a))

    e1 = GB.Line(bottom_left, bottom_right)
    e2 = GB.Line(bottom_right, top_right)
    e3 = GB.Line(top_right, top_left)
    e4 = GB.Line(top_left, bottom_left)

    return (e1, e2, e3, e4)
end

function get_vertices(a::GB.Rect, pos, axes)
    half_widths = maximum(a)
    half_width = half_widths[1]
    half_height = half_widths[2]

    x_cap = get_x_cap(axes)
    y_cap = get_y_cap(axes)

    bottom_left = pos .- half_width .* x_cap .- half_height .* y_cap
    bottom_right = pos .+ half_width .* x_cap .- half_height .* y_cap
    top_right = pos .+ half_width .* x_cap .+ half_height .* y_cap
    top_left = pos .- half_width .* x_cap .+ half_height .* y_cap

    return (bottom_left, bottom_right, top_right, top_left)
end

function get_vertices(a::GB.Rect)
    half_widths = maximum(a)
    half_width = half_widths[1]
    half_height = half_widths[2]

    VecType = typeof(half_widths)

    bottom_left = VecType(-half_width, -half_height)
    bottom_right = VecType(half_width, -half_height)
    top_right = VecType(half_width, half_height)
    top_left = VecType(-half_width, half_height)

    return (bottom_left, bottom_right, top_right, top_left)
end

function get_normals(a::GB.Rect)
    T = eltype(a.origin)
    axes = Axes{T}()
    x_cap = get_x_cap(axes)
    y_cap = get_y_cap(axes)
    return (-y_cap, x_cap, y_cap, -x_cap)
end

function get_normals(a::GB.Rect, pos, axes)
    x_cap = get_x_cap(axes)
    y_cap = get_y_cap(axes)
    return (-y_cap, x_cap, y_cap, -x_cap)
end
