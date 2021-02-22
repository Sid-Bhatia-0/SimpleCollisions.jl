#####
# Vec
#####

get_x(vec::SA.SVector) = vec[1]
get_y(vec::SA.SVector) = vec[2]

#####
# StdShape
#####

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

get_head(line::StdLine{T}) where {T} = SA.SVector(get_half_length(line), zero(T))
get_tail(line::StdLine{T}) where {T} = SA.SVector(-get_half_length(line), zero(T))
get_vertices(line::StdLine) = (get_head(line), get_tail(line))

function get_vertices(line::StdLine{T}, pos::SA.SVector{2, T}, axes::Axes{T}) where {T}
    half_length = get_half_length(line)
    x_cap = get_x_cap(axes)
    tail = pos .+ half_length .* -x_cap
    head = pos .+ half_length .* x_cap
    return (tail, head)
end

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
get_width(rect::StdRect) = 2 * get_half_width(rect)
get_half_height(rect::StdRect) = rect.half_height
get_height(rect::StdRect) = 2 * get_half_height(rect)

get_bottom_left(rect::StdRect) = SA.SVector(-get_half_width(rect), -get_half_height(rect))
get_bottom_right(rect::StdRect) = SA.SVector(get_half_width(rect), -get_half_height(rect))
get_top_right(rect::StdRect) = SA.SVector(get_half_width(rect), get_half_height(rect))
get_top_left(rect::StdRect) = SA.SVector(-get_half_width(rect), get_half_height(rect))

get_vertices(rect::StdRect{T}) where {T} = (get_bottom_left(rect), get_bottom_right(rect), get_top_right(rect), get_top_left(rect))

get_bottom_left(rect::StdRect{T}, pos::SA.SVector{2, T}) where {T} = pos .+ SA.SVector(-get_half_width(rect), -get_half_height(rect))
get_bottom_right(rect::StdRect{T}, pos::SA.SVector{2, T}) where {T} = pos .+ SA.SVector(get_half_width(rect), -get_half_height(rect))
get_top_right(rect::StdRect{T}, pos::SA.SVector{2, T}) where {T} = pos .+ SA.SVector(get_half_width(rect), get_half_height(rect))
get_top_left(rect::StdRect{T}, pos::SA.SVector{2, T}) where {T} = pos .+ SA.SVector(-get_half_width(rect), get_half_height(rect))

get_vertices(rect::StdRect{T}, pos::SA.SVector{2, T}) where {T} = (get_bottom_left(rect, pos), get_bottom_right(rect, pos), get_top_right(rect, pos), get_top_left(rect, pos))

function get_vertices(rect::StdRect{T}, pos::SA.SVector{2, T}, axes::Axes{T}) where {T}
    half_width = get_half_width(rect)
    half_height = get_half_height(rect)

    x_cap = get_x_cap(axes)
    y_cap = get_y_cap(axes)

    bottom_left = pos .- half_width .* x_cap .- half_height .* y_cap
    bottom_right = pos .+ half_width .* x_cap .- half_height .* y_cap
    top_right = pos .+ half_width .* x_cap .+ half_height .* y_cap
    top_left = pos .- half_width .* x_cap .+ half_height .* y_cap

    return (bottom_left, bottom_right, top_right, top_left)
end

get_normals(rect::StdRect{T}) where {T} = get_normals(rect, Axes{T}())

function get_normals(rect::StdRect{T}, axes::Axes{T}) where {T}
    x_cap = get_x_cap(axes)
    y_cap = get_y_cap(axes)
    return (-y_cap, x_cap, y_cap, -x_cap)
end

get_area(rect::StdRect{T}) where {T} = convert(T, 4 * get_half_width(rect) * get_half_height(rect))

function get_area(vertices)
    area = zero(eltype(vertices[1]))

    for i in 1:length(vertices) - 1
        v1 = augmented_vertices[i]
        v2 = augmented_vertices[i + 1]
        area = area + v1[1] * v2[2] - v2[1] * v1[2]
    end

    last_vertex = vertices[end]
    first_vertex = vertices[1]

    area = (area + last_vertex[1] * first_vertex[2] - first_vertex[1] * last_vertex[2]) / 2

    return area
end
