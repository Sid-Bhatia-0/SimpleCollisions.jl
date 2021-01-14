#####
# MassData
#####

struct MassData{T}
    mass::T
    inv_mass::T
end

function MassData(density::T, shape::GB.AbstractGeometry) where {T<:AbstractFloat}
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

function InertiaData(density::T, shape::GB.AbstractGeometry) where {T<:AbstractFloat}
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

mutable struct LinearMotionData{N, T}
    position::GB.Vec{N, T}
    velocity::GB.Vec{N, T}
end

LinearMotionData{N, T}() where {N, T} = LinearMotionData(zero(GB.Vec{N, T}), zero(GB.Vec{N, T}))

get_velocity(linear_motion_data::LinearMotionData) = linear_motion_data.velocity
set_velocity!(linear_motion_data::LinearMotionData, velocity) = linear_motion_data.velocity = velocity

get_position(linear_motion_data::LinearMotionData) = linear_motion_data.position
set_position!(linear_motion_data::LinearMotionData, position) = linear_motion_data.position = position

const LinearMotionData2{T} = LinearMotionData{2, T}

#####
# AngularMotionData2
#####

mutable struct AngularMotionData2{T}
    angle::T
    angular_velocity::T
end

AngularMotionData2{T}() where {T} = AngularMotionData2(zero(T), zero(T))

get_angle(angular_motion_data::AngularMotionData2) = angular_motion_data.angle
set_angle!(angular_motion_data::AngularMotionData2, angle) = angular_motion_data.angle = angle

get_angular_velocity(angular_motion_data::AngularMotionData2) = angular_motion_data.angular_velocity
set_angular_velocity!(angular_motion_data::AngularMotionData2, angular_velocity) = angular_motion_data.angular_velocity = angular_velocity

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
    linear_motion_data::LinearMotionData2{T}
    angular_motion_data::AngularMotionData2{T}
    force::GB.Vec2{T}
    torque::T
end

function RigidBody{T}() where {T<:AbstractFloat}
    shape = GB.Circle(zero(GB.Point2{T}), one(T))
    material_data = MaterialData{T}()
    mass_data = MassData(material_data.density, shape)
    inertia_data = InertiaData(material_data.density, shape)
    linear_motion_data = LinearMotionData2{T}()
    angular_motion_data = AngularMotionData2{T}()
    force = zero(GB.Vec2{T})
    torque = zero(T)

    return RigidBody{T, typeof(shape)}(shape, material_data, mass_data, inertia_data, linear_motion_data, angular_motion_data, force, torque)
end

MacroTools.@forward RigidBody.linear_motion_data get_position, set_position!, get_velocity, set_velocity!
MacroTools.@forward RigidBody.angular_motion_data get_angle, set_angle!, get_angular_velocity, set_angular_velocity!
MacroTools.@forward RigidBody.mass_data get_mass, get_inv_mass
MacroTools.@forward RigidBody.inertia_data get_inertia, get_inv_inertia
MacroTools.@forward RigidBody.material_data get_density, get_restitution

get_shape(body::RigidBody) = body.shape
set_shape!(body::RigidBody, shape) = body.shape = shape
get_shape(shape::GB.HyperSphere, position) = typeof(shape)(GB.Point(position), shape.r)
get_shape(shape::GB.HyperRectangle, position) = typeof(shape)(position..., shape.widths)

get_force(body::RigidBody) = body.force
set_force!(body::RigidBody, force) = body.force = force
get_torque(body::RigidBody) = body.torque
set_torque!(body::RigidBody, torque) = body.torque = torque

function update!(body::RigidBody, dt)
    inv_mass = get_inv_mass(body)
    force = get_force(body)
    velocity = get_velocity(body)
    new_velocity = velocity + inv_mass * force * dt
    set_velocity!(body, new_velocity)

    position = get_position(body)
    new_position = position + new_velocity * dt
    set_position!(body, new_position)

    shape = get_shape(body)
    new_shape = get_shape(shape, new_position)
    set_shape!(body, new_shape)

    return body
end
