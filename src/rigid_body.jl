#####
# MaterialData
#####

struct MaterialData{T<:AbstractFloat}
    density::T
    restitution::T
    static_friction_coeff::T
    dynamic_friction_coeff::T
end

MaterialData{T}() where {T<:AbstractFloat} = MaterialData(one(T), one(T), convert(T, 0.6), convert(T, 0.4))

get_density(material_data::MaterialData) = material_data.density
get_restitution(material_data::MaterialData) = material_data.restitution
get_static_friction_coeff(material_data::MaterialData) = material_data.static_friction_coeff
get_dynamic_friction_coeff(material_data::MaterialData) = material_data.dynamic_friction_coeff

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
# IntertiaData
#####

struct InertiaData{T}
    inertia::T
    inv_inertia::T
end

InertiaData{T}() where {T} = InertiaData(one(T), one(T))

get_inertia(inertia_data::InertiaData) = inertia_data.inertia
get_inv_inertia(inertia_data::InertiaData) = inertia_data.inv_inertia

function InertiaData(density::T, shape::GB.AbstractGeometry) where {T<:AbstractFloat}
    inertia = get_inertia(density, shape)
    inv_inertia = 1 / inertia
    InertiaData(convert(T, inertia), convert(T, inv_inertia))
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
    inertia_data::InertiaData{T}
    angle_accumulator::Accumulator{T}
    angular_velocity_accumulator::Accumulator{T}
    torque_accumulator::Accumulator{T}
    direction::GB.Vec2{T}
end

function RigidBody{T}() where {T<:AbstractFloat}
    shape = GB.HyperSphere(zero(GB.Point2{T}), one(T))
    material_data = MaterialData{T}()
    mass_data = MassData(material_data.density, shape)
    position_accumulator = Accumulator(zero(GB.Vec2{T}), zero(GB.Vec2{T}))
    velocity_accumulator = Accumulator(zero(GB.Vec2{T}), zero(GB.Vec2{T}))
    force_accumulator = Accumulator(zero(GB.Vec2{T}), zero(GB.Vec2{T}))
    inertia_data = InertiaData(material_data.density, shape)
    angle = zero(T)
    angle_accumulator = Accumulator(angle)
    angular_velocity_accumulator = Accumulator(zero(T))
    torque_accumulator = Accumulator(zero(T))
    direction = GB.Vec2{T}(cos(angle), sin(angle))

    return RigidBody(shape, material_data, mass_data, position_accumulator, velocity_accumulator, force_accumulator, inertia_data, angle_accumulator, angular_velocity_accumulator, torque_accumulator, direction)
end

get_shape(body::RigidBody) = body.shape
set_shape!(body::RigidBody, shape) = body.shape = shape

MacroTools.@forward RigidBody.material_data get_density, get_restitution, get_static_friction_coeff, get_dynamic_friction_coeff

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

MacroTools.@forward RigidBody.inertia_data get_inertia, get_inv_inertia

get_angle(body::RigidBody) = get_value(body.angle_accumulator)
add_angle_change!(body::RigidBody, change) = add_change!(body.angle_accumulator, change)
apply_angle_change!(body::RigidBody) = apply_change!(body.angle_accumulator)

get_angular_velocity(body::RigidBody) = get_value(body.angular_velocity_accumulator)
add_angular_velocity_change!(body::RigidBody, change) = add_change!(body.angular_velocity_accumulator, change)
apply_angular_velocity_change!(body::RigidBody) = apply_change!(body.angular_velocity_accumulator)

get_torque(body::RigidBody) = get_value(body.torque_accumulator)
add_torque_change!(body::RigidBody, change) = add_change!(body.torque_accumulator, change)
apply_torque_change!(body::RigidBody) = apply_change!(body.torque_accumulator)

get_direction(body::RigidBody) = body.direction

@pretty_print RigidBody
