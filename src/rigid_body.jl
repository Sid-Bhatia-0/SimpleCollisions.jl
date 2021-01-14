#####
# MassData
#####

struct MassData{T}
    mass::T
    inv_mass::T
end

function MassData(density::T, shape::GB.GeometryPrimitive{2}) where {T<:AbstractFloat}
    mass = get_mass(density, shape)
    inv_mass = 1 / mass
    MassData(convert(T, mass), convert(T, inv_mass))
end

MassData{T}() where {T} = MassData(one(T), one(T))

get_mass(mass_data::MassData) = mass_data.mass
get_inv_mass(mass_data::MassData) = mass_data.inv_mass

#####
# IntertiaData
#####

struct InertiaData{T}
    inertia::T
    inv_inertia::T
end

function InertiaData(density::T, shape::GB.GeometryPrimitive{2}) where {T<:AbstractFloat}
    inertia = get_inertia(density, shape)
    inv_inertia = 1 / inertia
    InertiaData(convert(T, inertia), convert(T, inv_inertia))
end

InertiaData{T}() where {T} = InertiaData(one(T), one(T))

get_inertia(inertia_data::InertiaData) = inertia_data.inertia
get_inv_inertia(inertia_data::InertiaData) = inertia_data.inv_inertia

#####
# LinearMotionData
#####

mutable struct LinearMotionData{T}
    position::GB.Vec2{T}
    velocity::GB.Vec2{T}
end

LinearMotionData{T}() where {T} = LinearMotionData(zero(GB.Vec2{T}), zero(GB.Vec2{T}))

get_velocity(linear_motion_data::LinearMotionData) = linear_motion_data.velocity
set_velocity!(linear_motion_data::LinearMotionData, velocity) = linear_motion_data.velocity = velocity

get_position(linear_motion_data::LinearMotionData) = linear_motion_data.position
set_position!(linear_motion_data::LinearMotionData, position) = linear_motion_data.position = position

#####
# AngularMotionData
#####

mutable struct AngularMotionData{T}
    angle::T
    angular_velocity::T
end

AngularMotionData{T}() where {T} = AngularMotionData(zero(T), zero(T))

get_angle(angular_motion_data::AngularMotionData) = angular_motion_data.angle
set_angle!(angular_motion_data::AngularMotionData, angle) = angular_motion_data.angle = angle

get_angular_velocity(angular_motion_data::AngularMotionData) = angular_motion_data.angular_velocity
set_angular_velocity!(angular_motion_data::AngularMotionData, angular_velocity) = angular_motion_data.angular_velocity = angular_velocity

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
# RigidBody
#####

mutable struct RigidBody{T<:AbstractFloat, S<:GB.AbstractGeometry{2, T}}
    shape::S
    material_data::MaterialData{T}
    mass_data::MassData{T}
    inertia_data::InertiaData{T}
    linear_motion_data::LinearMotionData{T}
    angular_motion_data::AngularMotionData{T}
    force::GB.Vec2{T}
    torque::T
end

function RigidBody{T}() where {T<:AbstractFloat}
    shape = GB.Circle(zero(GB.Point2{T}), one(T))
    material_data = MaterialData{T}()
    mass_data = MassData(material_data.density, shape)
    inertia_data = InertiaData(material_data.density, shape)
    linear_motion_data = LinearMotionData{T}()
    angular_motion_data = AngularMotionData{T}()
    force = zero(GB.Vec2{T})
    torque = zero(T)

    return RigidBody(shape, material_data, mass_data, inertia_data, linear_motion_data, angular_motion_data, force, torque)
end

MacroTools.@forward RigidBody.linear_motion_data get_position, set_position!, get_velocity, set_velocity!
MacroTools.@forward RigidBody.angular_motion_data get_angle, set_angle!, get_angular_velocity, set_angular_velocity!
MacroTools.@forward RigidBody.mass_data get_mass, get_inv_mass
MacroTools.@forward RigidBody.inertia_data get_inertia, get_inv_inertia
MacroTools.@forward RigidBody.material_data get_density, get_restitution
