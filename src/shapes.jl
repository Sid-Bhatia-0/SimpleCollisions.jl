abstract type AbstractStdShape{T} end

#####
# StdPoint
#####

struct StdPoint{T} <: AbstractStdShape{T} end

#####
# StdLine
#####

struct StdLine{T} <: AbstractStdShape{T}
    half_length::T
end

get_half_length(line::StdLine) = line.half_length

# head, tail, and vertices for StdLine
get_head(line::StdLine{T}) where {T} = Vector2D(get_half_length(line), zero(T))
get_tail(line::StdLine{T}) where {T} = Vector2D(-get_half_length(line), zero(T))
get_vertices(line::StdLine) = (get_tail(line), get_head(line))

# head, tail, and vertices for StdLine at arbitrary position
get_head(line::StdLine{T}, pos::Vector2D{T}) where {T} = pos + get_head(line)
get_tail(line::StdLine{T}, pos::Vector2D{T}) where {T} = pos + get_tail(line)
get_vertices(line::StdLine{T}, pos::Vector2D{T}) where {T} = (get_tail(line, pos), get_head(line, pos))

# head, tail, and vertices for StdLine at arbitrary position and orientation
get_head(line::StdLine{T}, pos::Vector2D{T}, dir::Vector2D{T}) where {T} = pos + get_half_length(line) * dir
get_tail(line::StdLine{T}, pos::Vector2D{T}, dir::Vector2D{T}) where {T} = pos - get_half_length(line) * dir
get_vertices(line::StdLine{T}, pos::Vector2D{T}, dir::Vector2D{T}) where {T} = (get_tail(line, pos, dir), get_head(line, pos, dir))

#####
# StdCircle
#####

struct StdCircle{T} <: AbstractStdShape{T}
    radius::T
end

get_radius(circle::StdCircle) = circle.radius

function get_area(circle::StdCircle{T}) where {T}
    radius = get_radius(circle)
    return convert(T, pi * radius * radius)
end

#####
# StdRect
#####

struct StdRect{T} <: AbstractStdShape{T}
    half_width::T
    half_height::T
end

get_half_width(rect::StdRect) = rect.half_width
get_half_height(rect::StdRect) = rect.half_height

get_width(rect::StdRect) = 2 * get_half_width(rect)
get_height(rect::StdRect) = 2 * get_half_height(rect)

# bottom_left, bottom_right, top_left, top_right, and vertices for StdRect
get_bottom_left(rect::StdRect) = Vector2D(-get_half_width(rect), -get_half_height(rect))
get_bottom_right(rect::StdRect) = Vector2D(get_half_width(rect), -get_half_height(rect))
get_top_right(rect::StdRect) = Vector2D(get_half_width(rect), get_half_height(rect))
get_top_left(rect::StdRect) = Vector2D(-get_half_width(rect), get_half_height(rect))
get_vertices(rect::StdRect{T}) where {T} = (get_bottom_left(rect), get_bottom_right(rect), get_top_right(rect), get_top_left(rect))

# bottom_left, bottom_right, top_left, top_right, and vertices for StdRect at arbitrary position
get_bottom_left(rect::StdRect{T}, pos::Vector2D{T}) where {T} = pos + Vector2D(-get_half_width(rect), -get_half_height(rect))
get_bottom_right(rect::StdRect{T}, pos::Vector2D{T}) where {T} = pos + Vector2D(get_half_width(rect), -get_half_height(rect))
get_top_right(rect::StdRect{T}, pos::Vector2D{T}) where {T} = pos + Vector2D(get_half_width(rect), get_half_height(rect))
get_top_left(rect::StdRect{T}, pos::Vector2D{T}) where {T} = pos + Vector2D(-get_half_width(rect), get_half_height(rect))
get_vertices(rect::StdRect{T}, pos::Vector2D{T}) where {T} = (get_bottom_left(rect, pos), get_bottom_right(rect, pos), get_top_right(rect, pos), get_top_left(rect, pos))

# bottom_left, bottom_right, top_left, top_right, and vertices for StdRect at arbitrary position and orientation
get_bottom_left(rect::StdRect{T}, pos::Vector2D{T}, dir::Vector2D{T}) where {T} = pos - get_half_width(rect) * dir - get_half_height(rect) * rotate_plus_90(dir)
get_bottom_right(rect::StdRect{T}, pos::Vector2D{T}, dir::Vector2D{T}) where {T} = pos + get_half_width(rect) * dir - get_half_height(rect) * rotate_plus_90(dir)
get_top_right(rect::StdRect{T}, pos::Vector2D{T}, dir::Vector2D{T}) where {T} = pos + get_half_width(rect) * dir + get_half_height(rect) * rotate_plus_90(dir)
get_top_left(rect::StdRect{T}, pos::Vector2D{T}, dir::Vector2D{T}) where {T} = pos - get_half_width(rect) * dir + get_half_height(rect) * rotate_plus_90(dir)
get_vertices(rect::StdRect{T}, pos::Vector2D{T}, dir::Vector2D{T}) where {T} = (get_bottom_left(rect, pos, dir), get_bottom_right(rect, pos, dir), get_top_right(rect, pos, dir), get_top_left(rect, pos, dir))

get_area(rect::StdRect{T}) where {T} = convert(T, 4 * get_half_width(rect) * get_half_height(rect))

get_normals(rect::StdRect{T}) where {T} = get_normals(rect, Vector2D(one(T), zero(T)))

function get_normals(rect::StdRect{T}, dir::Vector2D{T}) where {T}
    x_cap = dir
    y_cap = rotate_plus_90(dir)
    return (-y_cap, x_cap, y_cap, -x_cap)
end
