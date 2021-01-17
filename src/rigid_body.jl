#####
# MaterialData
#####

struct MaterialData{T}
    density::T
    restitution::T
end

MaterialData{T}() where {T} = MaterialData(one(T), one(T))

get_density(material_data::MaterialData) = material_data.density
get_restitution(material_data::MaterialData) = material_data.restitution

#####
# MassData
#####

struct MassData{T}
    mass::T
    inv_mass::T
end

MassData{T}() where {T} = MassData(one(T), one(T))

get_mass(mass_data::MassData) = mass_data.mass
get_inv_mass(mass_data::MassData) = mass_data.inv_mass

function MassData(density::T, shape::GB.AbstractGeometry) where {T<:AbstractFloat}
    mass = get_mass(density, shape)
    inv_mass = 1 / mass
    MassData(convert(T, mass), convert(T, inv_mass))
end

#####
# Accumulator
#####

mutable struct Accumulator{T}
    value::T
    change::T
end

Accumulator{T}() where {T} = Accumulator(one(T), one(T))

get_value(accumulator::Accumulator) = accumulator.value
set_value!(accumulator::Accumulator, value) = accumulator.value = value

get_change(accumulator::Accumulator) = accumulator.change
set_change!(accumulator::Accumulator, change) = accumulator.change = change

add_change!(accumulator::Accumulator, change) = set_change!(accumulator, get_change(accumulator) .+ change)

function apply_change!(accumulator::Accumulator{T}) where {T}
    set_value!(accumulator, get_value(accumulator) .+ get_change(accumulator))
    set_change!(accumulator, zero(T))
end

#####
# RigidBody
#####

mutable struct RigidBody{T<:AbstractFloat, S<:GB.AbstractGeometry{2, T}}
    shape::S
    material_data::MaterialData{T}
    mass_data::MassData{T}
    position_accumulator::Accumulator{GB.Vec2{T}}
    velocity_accumulator::Accumulator{GB.Vec2{T}}
    force_accumulator::Accumulator{GB.Vec2{T}}
end

function RigidBody{T}() where {T<:AbstractFloat}
    shape = GB.HyperSphere(zero(GB.Point2{T}), one(T))
    material_data = MaterialData{T}()
    mass_data = MassData(material_data.density, shape)
    position_accumulator = Accumulator(zero(GB.Vec2{T}), zero(GB.Vec2{T}))
    velocity_accumulator = Accumulator(zero(GB.Vec2{T}), zero(GB.Vec2{T}))
    force_accumulator = Accumulator(zero(GB.Vec2{T}), zero(GB.Vec2{T}))

    return RigidBody(shape, material_data, mass_data, position_accumulator, velocity_accumulator, force_accumulator)
end

get_shape(body::RigidBody) = body.shape
set_shape!(body::RigidBody, shape) = body.shape = shape

MacroTools.@forward RigidBody.material_data get_density, get_restitution

MacroTools.@forward RigidBody.mass_data get_mass, get_inv_mass

get_position(body::RigidBody) = get_value(body.position_accumulator)
add_position_change!(body::RigidBody, change) = add_change!(body.position_accumulator, change)
apply_position_change!(body::RigidBody) = apply_change!(body.position_accumulator)

get_velocity(body::RigidBody) = get_value(body.velocity_accumulator)
add_velocity_change!(body::RigidBody, change) = add_change!(body.velocity_accumulator, change)
apply_velocity_change!(body::RigidBody) = apply_change!(body.velocity_accumulator)

get_force(body::RigidBody) = get_value(body.force_accumulator)
add_force_change!(body::RigidBody, change) = add_change!(body.force_accumulator, change)
apply_force_change!(body::RigidBody) = apply_change!(body.force_accumulator)

translated(shape::GB.HyperSphere, position) = typeof(shape)(GB.Point(position), shape.r)
translated(shape::GB.HyperRectangle, position) = typeof(shape)(position..., shape.widths...)

@pretty_print RigidBody
