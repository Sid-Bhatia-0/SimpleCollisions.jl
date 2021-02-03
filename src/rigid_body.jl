#####
# MaterialData
#####

struct MaterialData{T<:AbstractFloat}
    density::T
    restitution::T
    static_friction_coeff::T
    dynamic_friction_coeff::T
end

get_density(material_data::MaterialData) = material_data.density
get_restitution(material_data::MaterialData) = material_data.restitution
get_static_friction_coeff(material_data::MaterialData) = material_data.static_friction_coeff
get_dynamic_friction_coeff(material_data::MaterialData) = material_data.dynamic_friction_coeff

MaterialData{T}() where {T<:AbstractFloat} = MaterialData(one(T), one(T), convert(T, 0.6), convert(T, 0.4))

#####
# MassData
#####

struct MassData{T<:AbstractFloat}
    mass::T
    inv_mass::T
end

get_mass(mass_data::MassData) = mass_data.mass
get_inv_mass(mass_data::MassData) = mass_data.inv_mass

MassData{T}() where {T<:AbstractFloat} = MassData(one(T), one(T))

get_mass(density, shape) = density * get_area(shape)

function MassData(density::T, shape::GB.AbstractGeometry) where {T<:AbstractFloat}
    mass = get_mass(density, shape)
    inv_mass = 1 / mass
    return MassData(convert(T, mass), convert(T, inv_mass))
end

#####
# IntertiaData
#####

struct InertiaData{T<:AbstractFloat}
    inertia::T
    inv_inertia::T
end

get_inertia(inertia_data::InertiaData) = inertia_data.inertia
get_inv_inertia(inertia_data::InertiaData) = inertia_data.inv_inertia

InertiaData{T}() where {T<:AbstractFloat} = InertiaData(one(T), one(T))

get_inertia(density, shape::GB.Circle) = get_mass(density, shape) * shape.r * shape.r / 2
get_inertia(density, shape::GB.Rect2D) = get_mass(density, shape) * LA.dot(shape.widths, shape.widths) / 12

function InertiaData(density::T, shape::GB.AbstractGeometry) where {T<:AbstractFloat}
    inertia = get_inertia(density, shape)
    inv_inertia = 1 / inertia
    InertiaData(convert(T, inertia), convert(T, inv_inertia))
end

#####
# Axes
#####

mutable struct Axes{T}
    x_cap::GB.Vec2{T}
    y_cap::GB.Vec2{T}
end

get_x_cap(axes::Axes) = axes.x_cap
set_x_cap!(axes::Axes, x_cap) = axes.x_cap = x_cap
get_y_cap(axes::Axes) = axes.y_cap
set_y_cap!(axes::Axes, y_cap) = axes.y_cap = y_cap

function Axes{T}() where {T}
    x_cap = GB.unit(GB.Vec2{T}, 1)
    y_cap = GB.unit(GB.Vec2{T}, 2)
    return Axes{T}(x_cap, y_cap)
end

function Axes(angle::T) where {T}
    x_cap = GB.Vec2{T}(cos(angle), sin(angle))
    y_cap = rotate_90(x_cap)
    return Axes{T}(x_cap, y_cap)
end

function rotate(vec::GB.Vec2, axes::Axes)
    x_cap = get_x_cap(axes)
    y_cap = get_y_cap(axes)
    v1 = vec[1]
    v2 = vec[2]
    return typeof(vec)(x_cap[1] * v1 + y_cap[1] * v2, x_cap[2] * v1 + y_cap[2] * v2)
end

rotate_90(vec::GB.Vec2) = typeof(vec)(-vec[2], vec[1])
rotate_180(vec::GB.Vec2) = -vec
rotate_minus_90(vec::GB.Vec2) = typeof(vec)(vec[2], -vec[1])

function rotate_90(axes::Axes)
    x_cap = get_x_cap(axes)
    y_cap = get_y_cap(axes)
    x_cap_90 = y_cap
    y_cap_90 = -x_cap
    return Axes(x_cap_90, y_cap_90)
end

function rotate_180(axes::Axes)
    x_cap = get_x_cap(axes)
    y_cap = get_y_cap(axes)
    return Axes(-x_cap, -y_cap)
end

function rotate_minus_90(axes::Axes)
    x_cap = get_x_cap(axes)
    y_cap = get_y_cap(axes)
    x_cap_minus_90 = -y_cap
    y_cap_minus_90 = x_cap
    return Axes(x_cap_minus_90, y_cap_minus_90)
end

get_relative_direction(d1::GB.Vec2, d2::GB.Vec2) = typeof(d1)(d2[1] * d1[1] + d2[2] * d1[2], d2[2] * d1[1] - d2[1] * d1[2])
invert_relative_direction(d::GB.Vec) = typeof(d)(d[1], -d[2])

function get_relative_axes(axes1::Axes, axes2::Axes)
    x_cap_21 = get_relative_direction(get_x_cap(axes1), get_x_cap(axes2))
    y_cap_21 = rotate_90(x_cap_21)
    return Axes(x_cap_21, y_cap_21)
end

function invert_relative_axes(axes::Axes)
    x_cap_21 = get_x_cap(axes)
    y_cap_21 = get_y_cap(axes)
    x_cap_12 = invert_relative_direction(x_cap_21)
    y_cap_12 = rotate_90(x_cap_12)
    return Axes(x_cap_12, y_cap_12)
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

struct RigidBody{T<:AbstractFloat, S<:GB.AbstractGeometry{2, T}}
    shape::S
    material_data::MaterialData{T}

    mass_data::MassData{T}
    position_accumulator::Accumulator{GB.Vec2{T}}
    velocity_accumulator::Accumulator{GB.Vec2{T}}
    force_accumulator::Accumulator{GB.Vec2{T}}

    inertia_data::InertiaData{T}
    angle_accumulator::Accumulator{T}
    axes::Axes{T}
    angular_velocity_accumulator::Accumulator{T}
    torque_accumulator::Accumulator{T}
end

function RigidBody{T}() where {T<:AbstractFloat}
    PointType = GB.Point2{T}
    VecType = GB.Vec2{T}

    shape = GB.HyperSphere(zero(PointType), one(T))
    material_data = MaterialData{T}()

    mass_data = MassData(material_data.density, shape)
    position_accumulator = Accumulator(zero(VecType), zero(VecType))
    velocity_accumulator = Accumulator(zero(VecType), zero(VecType))
    force_accumulator = Accumulator(zero(VecType), zero(VecType))

    inertia_data = InertiaData(material_data.density, shape)
    angle = zero(T)
    angle_accumulator = Accumulator(angle, zero(T))
    axes = Axes(angle)
    angular_velocity_accumulator = Accumulator(zero(T), zero(T))
    torque_accumulator = Accumulator(zero(T), zero(T))

    return RigidBody(shape,
                     material_data,
                     mass_data,
                     position_accumulator,
                     velocity_accumulator,
                     force_accumulator,
                     inertia_data,
                     angle_accumulator,
                     axes,
                     angular_velocity_accumulator,
                     torque_accumulator)
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

get_axes(body::RigidBody) = body.axes
MacroTools.@forward RigidBody.axes get_x_cap, set_x_cap!, get_y_cap, set_y_cap!

get_angular_velocity(body::RigidBody) = get_value(body.angular_velocity_accumulator)
add_angular_velocity_change!(body::RigidBody, change) = add_change!(body.angular_velocity_accumulator, change)
apply_angular_velocity_change!(body::RigidBody) = apply_change!(body.angular_velocity_accumulator)

get_torque(body::RigidBody) = get_value(body.torque_accumulator)
add_torque_change!(body::RigidBody, change) = add_change!(body.torque_accumulator, change)
apply_torque_change!(body::RigidBody) = apply_change!(body.torque_accumulator)

@pretty_print RigidBody
