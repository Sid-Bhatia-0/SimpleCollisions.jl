rotate_plus_90(vec::Vector2D) = typeof(vec)(-vec[2], vec[1])
rotate_minus_90(vec::Vector2D) = typeof(vec)(vec[2], -vec[1])
rotate_180(vec::Vector2D) = -vec

rotate(x, y, c, s) = Vector2D(c * x - s * y, s * x + c * y)
rotate(vec::Vector2D, theta) = rotate(vec[1], vec[2], cos(theta), sin(theta))
rotate(vec::Vector2D, dir::Vector2D) = rotate(vec[1], vec[2], dir[1], dir[2])

get_relative_direction(dir1::Vector2D, dir2::Vector2D) = typeof(dir1)(dir2[1] * dir1[1] + dir2[2] * dir1[2], dir2[2] * dir1[1] - dir2[1] * dir1[2])
invert_relative_direction(dir::Vector2D) = typeof(dir)(dir[1], -dir[2])

function invert(pos::Vector2D, dir::Vector2D)
    inv_dir = invert_relative_direction(dir)
    inv_pos = -rotate(pos, inv_dir)
    return inv_pos, inv_dir
end
