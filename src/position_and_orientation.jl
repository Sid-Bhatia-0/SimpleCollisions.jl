struct Vector2D{T}
    x::T
    y::T
end

Base.:-(vector::Vector2D) = Vector2D(-vector.x, -vector.y)
Base.:+(vector1::Vector2D, vector2::Vector2D) = Vector2D(vector1.x + vector2.x, vector1.y + vector2.y)
Base.:-(vector1::Vector2D, vector2::Vector2D) = Vector2D(vector1.x - vector2.x, vector1.y - vector2.y)
Base.:*(number::Number, vector::Vector2D) = Vector2D(number * vector.x, number * vector.y)
Base.:*(vector::Vector2D, number::Number) = Vector2D(number * vector.x, number * vector.y)
Base.:/(vector::Vector2D, number::Number) = Vector2D(vector.x / number, vector.y / number)

get_x(vec::Vector2D) = vec[1]
get_y(vec::Vector2D) = vec[2]

rotate_plus_90(vec::Vector2D) = typeof(vec)(-vec[2], vec[1])
rotate_minus_90(vec::Vector2D) = typeof(vec)(vec[2], -vec[1])
rotate_180(vec::Vector2D) = -vec

rotate(x::T, y::T, c::T, s::T) where {T} = Vector2D(c * x - s * y, s * x + c * y)
rotate(vec::Vector2D{T}, theta::T) where {T} = rotate(vec[1], vec[2], cos(theta), sin(theta))
rotate(vec::Vector2D{T}, dir::Vector2D{T}) where {T} = rotate(vec[1], vec[2], dir[1], dir[2])

get_relative_direction(dir1::Vector2D{T}, dir2::Vector2D{T}) where {T} = typeof(dir1)(dir2[1] * dir1[1] + dir2[2] * dir1[2], dir2[2] * dir1[1] - dir2[1] * dir1[2])
invert_relative_direction(dir::Vector2D) = typeof(dir)(dir[1], -dir[2])

function invert(pos::Vector2D{T}, dir::Vector2D{T}) where {T}
    inv_dir = invert_relative_direction(dir)
    inv_pos = -rotate(pos, inv_dir)
    return inv_pos, inv_dir
end

#####
##### linear algebra
#####

LA.norm(vector::Vector2D) = LA.norm(SA.SVector(vector.x, vector.y))

function LA.normalize(vector::Vector2D)
    svector = SA.SVector(vector.x, vector.y)
    svector_normed = LA.normalize(svector)
    vector_normed = Vector2D(svector_normed[1], svector_normed[2])
    return vector_normed
end
