#####
# MassData
#####

struct MassData{T}
    mass::T
    inv_mass::T
end

function MassData(density::T, shape::GB.GeometryPrimitive{2}) where {T<:Number}
    mass = get_mass(density, shape)
    inv_mass = 1 / mass
    MassData(mass, inv_mass)
end

MassData{T}() where {T} = MassData(one(T), one(T))

#####
# IntertiaData
#####

struct InertiaData{T}
    inertia::T
    inv_inertia::T
end

function InertiaData(density::T, shape::GB.GeometryPrimitive{2}) where {T<:Number}
    inertia = get_inertia(density, shape)
    inv_inertia = 1 / inertia
    InertiaData(inertia, inv_inertia)
end

InertiaData{T}() where {T} = InertiaData(one(T), one(T))

#####
# LinearMotionData
#####

struct LinearMotionData{T}
    position::GB.Vec2{T}
    velocity::GB.Vec2{T}
end

LinearMotionData{T}() where {T} = LinearMotionData(zero(GB.Vec2{T}), zero(GB.Vec2{T}))

#####
# AngularMotionData
#####

struct AngularMotionData{T}
    angle::T
    angular_velocity::T
end

AngularMotionData{T}() where {T} = AngularMotionData(zero(T), zero(T))

#####
# MaterialData
#####

struct MaterialData{T}
    density::T
    restitution::T
end

MaterialData{T}() where {T} = MaterialData(one(T), one(T))

#####
# RigidBody
#####

mutable struct RigidBody{T<:AbstractFloat}
    shape::GB.GeometryPrimitive{2, T}
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
