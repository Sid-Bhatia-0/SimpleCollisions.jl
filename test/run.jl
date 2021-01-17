import PhysicsEngine2D
import PhysicsEngine2D: PE2D
import GeometryBasics
const GB = GeometryBasics
import Makie

const T = Float32
const NUM_ITER = 500
const DT = 0.01
const FRAME_RATE = 1 / DT

shape1 = GB.HyperSphere(GB.Point2{T}(5.0f0, 1.0f0), one(T))
material_data1 = PE2D.MaterialData{T}()
mass_data1 = PE2D.MassData(material_data1.density, shape1)
position_accumulator1 = PE2D.Accumulator(GB.Vec2{T}(5.0f0, 1.0f0), zero(GB.Vec2{T}))
velocity_accumulator1 = PE2D.Accumulator(GB.Vec2{T}(0.0f0, 1.0f0), zero(GB.Vec2{T}))
force_accumulator1 = PE2D.Accumulator(zero(GB.Vec2{T}), zero(GB.Vec2{T}))
body1 = PE2D.RigidBody(shape1, material_data1, mass_data1, position_accumulator1, velocity_accumulator1, force_accumulator1)

shape2 = GB.HyperSphere(GB.Point2{T}(1.0f0, 5.0f0), one(T))
material_data2 = PE2D.MaterialData{T}()
mass_data2 = PE2D.MassData(material_data2.density, shape2)
position_accumulator2 = PE2D.Accumulator(GB.Vec2{T}(1.0f0, 5.0f0), zero(GB.Vec2{T}))
velocity_accumulator2 = PE2D.Accumulator(GB.Vec2{T}(1.0f0, 0.0f0), zero(GB.Vec2{T}))
force_accumulator2 = PE2D.Accumulator(zero(GB.Vec2{T}), zero(GB.Vec2{T}))
body2 = PE2D.RigidBody(shape2, material_data2, mass_data2, position_accumulator2, velocity_accumulator2, force_accumulator2)

bodies = [body1, body2]
world = PE2D.World(bodies)
display(world)

scene = PE2D.render!(world, num_iter = NUM_ITER, dt = DT, frame_rate = FRAME_RATE)
