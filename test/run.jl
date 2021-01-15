import PhysicsEngine2D
import PhysicsEngine2D: PE2D
import GeometryBasics
const GB = GeometryBasics
import Makie

const NUM_ITER = 100
const DT = 0.01
const FRAME_RATE = 24

body1 = PE2D.RigidBody{Float32}()
body2 = PE2D.RigidBody{Float32}()
PE2D.set_velocity!(body1, GB.Vec2(1.0f0, 0.0f0))
PE2D.set_velocity!(body2, GB.Vec2(0.0f0, 1.0f0))

bodies = [body1, body2]
world = PE2D.World(bodies)

scene = PE2D.render!(world, num_iter = NUM_ITER, dt = DT, frame_rate = FRAME_RATE)
