function is_colliding(a::GB.Rect2D, b::GB.Rect2D)
    if (a.origin[1] + a.widths[1] < b.origin[1]) || (b.origin[1] + b.widths[1] < a.origin[1])
        return false
    end
    if (a.origin[2] + a.widths[2] < b.origin[2]) || (b.origin[2] + b.widths[2] < a.origin[2])
        return false
    end
    return true
end

function is_colliding(a::GB.Circle, b::GB.Circle)
    ba = b.center .- a.center
    return LA.dot(ba, ba) <= (a.r + b.r) ^ 2
end
