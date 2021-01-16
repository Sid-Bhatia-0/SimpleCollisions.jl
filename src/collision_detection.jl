function is_colliding(a::GB.HyperRectangle{N}, b::GB.HyperRectangle{N}) where {N}
    for i in 1:N
        if (a.origin[i] + a.widths[i] < b.origin[i]) || (b.origin[i] + b.widths[i] < a.origin[i])
            return false
        end
    end
    return true
end

function is_colliding(a::GB.HyperSphere{N}, b::GB.Point{N}) where {N}
    ba = b .- a.center
    return LA.dot(ba, ba) <= a.r ^ 2
end

function is_colliding(a::GB.HyperSphere{N}, b::GB.HyperSphere{N}) where {N}
    ba = b.center .- a.center
    return LA.dot(ba, ba) <= (a.r + b.r) ^ 2
end
