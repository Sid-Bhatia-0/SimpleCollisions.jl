get_x(vec::SA.SVector) = vec[1]
get_y(vec::SA.SVector) = vec[2]

rotate_plus_90(vec::SA.SVector{2}) = typeof(vec)(-vec[2], vec[1])
rotate_minus_90(vec::SA.SVector{2}) = typeof(vec)(vec[2], -vec[1])
rotate_180(vec::SA.SVector{2}) = -vec

rotate(x::T, y::T, c::T, s::T) where {T} = SA.SVector(c * x - s * y, s * x + c * y)
rotate(vec::SA.SVector{2, T}, theta::T) where {T} = rotate(vec[1], vec[2], cos(theta), sin(theta))
rotate(vec::SA.SVector{2, T}, dir::SA.SVector{2, T}) where {T} = rotate(vec[1], vec[2], dir[1], dir[2])

get_relative_direction(dir1::SA.SVector{2, T}, dir2::SA.SVector{2, T}) where {T} = typeof(dir1)(dir2[1] * dir1[1] + dir2[2] * dir1[2], dir2[2] * dir1[1] - dir2[1] * dir1[2])
invert_relative_direction(dir::SA.SVector{2}) = typeof(dir)(dir[1], -dir[2])

function invert(pos::SA.SVector{2, T}, dir::SA.SVector{2, T}) where {T}
    inv_dir = invert_relative_direction(dir)
    inv_pos = -rotate(pos, inv_dir)
    return inv_pos, inv_dir
end
