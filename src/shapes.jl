abstract type AbstractShape end

#####
# StandardPoint
#####

struct StandardPoint{T} <: AbstractShape end

#####
# StandardLine
#####

struct StandardLine{T} <: AbstractShape
    half_length::T
end

get_half_length(line::StandardLine) = line.half_length

# head, tail, and vertices for StandardLine
get_head(line::StandardLine) = Vector2D(get_half_length(line), zero(get_half_length(line)))
get_tail(line::StandardLine) = Vector2D(-get_half_length(line), zero(get_half_length(line)))
get_vertices(line::StandardLine) = (get_tail(line), get_head(line))

# head, tail, and vertices for StandardLine at arbitrary position
get_head(line::StandardLine, pos) = pos + get_head(line)
get_tail(line::StandardLine, pos) = pos + get_tail(line)
get_vertices(line::StandardLine, pos) = (get_tail(line, pos), get_head(line, pos))

# head, tail, and vertices for StandardLine at arbitrary position and orientation
get_head(line::StandardLine, pos, dir) = pos + get_half_length(line) * dir
get_tail(line::StandardLine, pos, dir) = pos - get_half_length(line) * dir
get_vertices(line::StandardLine, pos, dir) = (get_tail(line, pos, dir), get_head(line, pos, dir))

#####
# StandardCircle
#####

struct StandardCircle{T} <: AbstractShape
    radius::T
end

get_radius(circle::StandardCircle) = circle.radius

function get_area(circle::StandardCircle)
    radius = get_radius(circle)
    return convert(typeof(radius), pi * radius * radius)
end

#####
# StandardRect
#####

struct StandardRect{T} <: AbstractShape
    half_width::T
    half_height::T
end

get_half_width(rect::StandardRect) = rect.half_width
get_half_height(rect::StandardRect) = rect.half_height

get_width(rect::StandardRect) = 2 * get_half_width(rect)
get_height(rect::StandardRect) = 2 * get_half_height(rect)

# bottom_left, bottom_right, top_left, top_right, and vertices for StandardRect
get_bottom_left(rect::StandardRect) = Vector2D(-get_half_width(rect), -get_half_height(rect))
get_bottom_right(rect::StandardRect) = Vector2D(get_half_width(rect), -get_half_height(rect))
get_top_right(rect::StandardRect) = Vector2D(get_half_width(rect), get_half_height(rect))
get_top_left(rect::StandardRect) = Vector2D(-get_half_width(rect), get_half_height(rect))
get_vertices(rect::StandardRect) = (get_bottom_left(rect), get_bottom_right(rect), get_top_right(rect), get_top_left(rect))

# bottom_left, bottom_right, top_left, top_right, and vertices for StandardRect at arbitrary position
get_bottom_left(rect::StandardRect, pos) = pos + typeof(pos)(-get_half_width(rect), -get_half_height(rect))
get_bottom_right(rect::StandardRect, pos) = pos + typeof(pos)(get_half_width(rect), -get_half_height(rect))
get_top_right(rect::StandardRect, pos) = pos + typeof(pos)(get_half_width(rect), get_half_height(rect))
get_top_left(rect::StandardRect, pos) = pos + typeof(pos)(-get_half_width(rect), get_half_height(rect))
get_vertices(rect::StandardRect, pos) = (get_bottom_left(rect, pos), get_bottom_right(rect, pos), get_top_right(rect, pos), get_top_left(rect, pos))

# bottom_left, bottom_right, top_left, top_right, and vertices for StandardRect at arbitrary position and orientation
get_bottom_left(rect::StandardRect, pos, dir) = pos - get_half_width(rect) * dir - get_half_height(rect) * rotate_plus_90(dir)
get_bottom_right(rect::StandardRect, pos, dir) = pos + get_half_width(rect) * dir - get_half_height(rect) * rotate_plus_90(dir)
get_top_right(rect::StandardRect, pos, dir) = pos + get_half_width(rect) * dir + get_half_height(rect) * rotate_plus_90(dir)
get_top_left(rect::StandardRect, pos, dir) = pos - get_half_width(rect) * dir + get_half_height(rect) * rotate_plus_90(dir)
get_vertices(rect::StandardRect, pos, dir) = (get_bottom_left(rect, pos, dir), get_bottom_right(rect, pos, dir), get_top_right(rect, pos, dir), get_top_left(rect, pos, dir))

get_area(rect::StandardRect) = convert(typeof(get_half_width(rect)), 4 * get_half_width(rect) * get_half_height(rect))

get_normals(rect::StandardRect) = get_normals(rect, Vector2D(one(typeof(get_half_width(rect))), zero(typeof(get_half_width(rect)))))

function get_normals(rect::StandardRect, dir)
    x_cap = dir
    y_cap = rotate_plus_90(dir)
    return (-y_cap, x_cap, y_cap, -x_cap)
end
