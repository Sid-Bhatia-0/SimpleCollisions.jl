#####
# AbstractBody
#####

abstract type AbstractBody{T} end

get_shape(body::AbstractBody) = body.shape

get_position(body::AbstractBody) = body.position
set_position!(body::AbstractBody{T}, position::SA.SVector{2, T}) where {T} = body.position = position

get_position_change(body::AbstractBody) = body.position_change
set_position_change!(body::AbstractBody{T}, position_change::SA.SVector{2, T}) where {T} = body.position_change = position_change

add_position_change!(body::AbstractBody{T}, position_change::SA.SVector{2, T}) where {T} = set_position_change!(body, get_position_change(body) .+ position_change)

function apply_position_change!(body)
    position = get_position(body)
    position_change = get_position_change(body)
    set_position!(body, position .+ position_change)
    set_position_change!(body, zero(position_change))
    return body
end

get_velocity(body::AbstractBody) = body.velocity
set_velocity!(body::AbstractBody{T}, velocity::SA.SVector{2, T}) where {T} = body.velocity = velocity

get_velocity_change(body::AbstractBody) = body.velocity_change
set_velocity_change!(body::AbstractBody{T}, velocity_change::SA.SVector{2, T}) where {T} = body.velocity_change = velocity_change

add_velocity_change!(body::AbstractBody{T}, velocity_change::SA.SVector{2, T}) where {T} = set_velocity_change!(body, get_velocity_change(body) .+ velocity_change)

function apply_velocity_change!(body)
    velocity = get_velocity(body)
    velocity_change = get_velocity_change(body)
    set_velocity!(body, velocity .+ velocity_change)
    set_velocity_change!(body, zero(velocity_change))
    return body
end

get_inv_mass(body::AbstractBody) = body.inv_mass
