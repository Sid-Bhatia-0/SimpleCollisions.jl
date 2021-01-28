import PhysicsEngine2D
import PhysicsEngine2D: PE2D
import GeometryBasics
const GB = GeometryBasics
import Makie

const T = Float32
const NUM_ITER = 500
const DT = 0.01
const FRAME_RATE = 1 / DT

shape1 = GB.HyperSphere(GB.Point2{T}(5, 1), one(T))
material_data1 = PE2D.MaterialData{T}()
mass_data1 = PE2D.MassData(material_data1.density, shape1)
position_accumulator1 = PE2D.Accumulator(GB.Vec2{T}(5, 1), zero(GB.Vec2{T}))
velocity_accumulator1 = PE2D.Accumulator(GB.Vec2{T}(0, 1), zero(GB.Vec2{T}))
force_accumulator1 = PE2D.Accumulator(zero(GB.Vec2{T}), zero(GB.Vec2{T}))
inertia_data1 = PE2D.InertiaData(material_data1.density, shape1)
angle_accumulator1 = PE2D.Accumulator(zero(T), zero(T))
angular_velocity_accumulator1 = PE2D.Accumulator(one(T), zero(T))
torque_accumulator1 = PE2D.Accumulator(zero(T), zero(T))
body1 = PE2D.RigidBody(shape1, material_data1, mass_data1, position_accumulator1, velocity_accumulator1, force_accumulator1, inertia_data1, angle_accumulator1, angular_velocity_accumulator1, torque_accumulator1)

shape2 = GB.HyperSphere(GB.Point2{T}(1, 5), one(T))
material_data2 = PE2D.MaterialData{T}()
mass_data2 = PE2D.MassData(material_data2.density, shape2)
position_accumulator2 = PE2D.Accumulator(GB.Vec2{T}(1, 5), zero(GB.Vec2{T}))
velocity_accumulator2 = PE2D.Accumulator(GB.Vec2{T}(1, 0), zero(GB.Vec2{T}))
force_accumulator2 = PE2D.Accumulator(zero(GB.Vec2{T}), zero(GB.Vec2{T}))
inertia_data2 = PE2D.InertiaData(material_data2.density, shape2)
angle_accumulator2 = PE2D.Accumulator(zero(T), zero(T))
angular_velocity_accumulator2 = PE2D.Accumulator(one(T), zero(T))
torque_accumulator2 = PE2D.Accumulator(zero(T), zero(T))
body2 = PE2D.RigidBody(shape2, material_data2, mass_data2, position_accumulator2, velocity_accumulator2, force_accumulator2, inertia_data2, angle_accumulator2, angular_velocity_accumulator2, torque_accumulator2)

bodies = [body1, body2]
world = PE2D.World(bodies)
display(world)

scene = PE2D.render!(world, num_iter = NUM_ITER, dt = DT, frame_rate = FRAME_RATE)
