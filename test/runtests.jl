import PhysicsEngine2D
import PhysicsEngine2D: PE2D
import GeometryBasics
const GB = GeometryBasics
using Test

function test_collision_list(collision_list)
    for (a, b, value) in collision_list
        @test PE2D.is_colliding(a, b) == value
        @test PE2D.is_colliding(b, a) == value
    end
end

function test_manifold_list(manifold_list)
    for (a, b, value) in manifold_list
        manifold = PE2D.Manifold(a, b)
        @test manifold.penetration ≈ value.penetration
        @test manifold.normal ≈ value.normal
    end
end

@testset "PhysicsEngine2D.jl" begin
    @testset "Area computation" begin
        @testset "Rect2D" begin
            a = GB.Rect(1, 2, 3, 4)
            @test GB.area(a) == 12
        end

        @testset "Circle" begin
            a = GB.HyperSphere(GB.Point(0, 0), 1)
            @test GB.area(a) ≈ π
        end
    end

    @testset "Collision detection" begin
        @testset "Point2 vs. Point2" begin
            collision_list = [(GB.Point(0, 0), GB.Point(0, 0), true),
                              (GB.Point(0, 0), GB.Point(1, 0), false)]
            test_collision_list(collision_list)
        end

        @testset "Circle vs. Point2" begin
            collision_list = [(GB.HyperSphere(GB.Point(0, 0), 1), GB.Point(0, 0), true),
                              (GB.HyperSphere(GB.Point(0, 0), 1), GB.Point(1, 0), true),
                              (GB.HyperSphere(GB.Point(0, 0), 1), GB.Point(2, 0), false)]
            test_collision_list(collision_list)
        end

        @testset "Circle vs. Circle" begin
            collision_list = [(GB.HyperSphere(GB.Point(0, 0), 1), GB.HyperSphere(GB.Point(0, 0), 1), true),
                              (GB.HyperSphere(GB.Point(0, 0), 1), GB.HyperSphere(GB.Point(1, 0), 1), true),
                              (GB.HyperSphere(GB.Point(0, 0), 1), GB.HyperSphere(GB.Point(2, 0), 1), true),
                              (GB.HyperSphere(GB.Point(0, 0), 1), GB.HyperSphere(GB.Point(3, 0), 1), false)]
            test_collision_list(collision_list)
        end

        @testset "Rect2D vs. Point2" begin
            collision_list = [(GB.Rect(1, 2, 5, 6), GB.Point(3, 3), true),
                              (GB.Rect(1, 2, 5, 6), GB.Point(1, 3), true),
                              (GB.Rect(1, 2, 5, 6), GB.Point(1, 2), true),
                              (GB.Rect(1, 2, 5, 6), GB.Point(1, 1), false)]
            test_collision_list(collision_list)
        end

        @testset "Rect2D vs. Circle" begin
            collision_list = [(GB.Rect(1, 2, 5, 6), GB.HyperSphere(GB.Point(3, 3), 1), true),
                              (GB.Rect(1, 2, 5, 6), GB.HyperSphere(GB.Point(4, 4), 1), true),
                              (GB.Rect(1, 2, 5, 6), GB.HyperSphere(GB.Point(4, 4), 10), true),
                              (GB.Rect(1, 2, 5, 6), GB.HyperSphere(GB.Point(0, 0), 3), true),
                              (GB.Rect(0, 1, 1, 1), GB.HyperSphere(GB.Point(0, 0), 1), true),
                              (GB.Rect(1, 2, 5, 6), GB.HyperSphere(GB.Point(0, 0), 1), false)]
            test_collision_list(collision_list)
        end

        @testset "Rect2D vs. Rect2D" begin
            collision_list = [(GB.Rect(1, 2, 3, 4), GB.Rect(3, 4, 5, 6), true),
                              (GB.Rect(1, 2, 3, 4), GB.Rect(4, 6, 1, 2), true),
                              (GB.Rect(1, 2, 3, 4), GB.Rect(0, 0, 6, 6), true),
                              (GB.Rect(1, 2, 3, 4), GB.Rect(4, 2, 1, 2), true),
                              (GB.Rect(1, 2, 3, 4), GB.Rect(5, 6, 7, 8), false)]
            test_collision_list(collision_list)
        end
    end

    @testset "Manifold generation" begin
        @testset "Circle vs. Circle" begin
            manifold_list = [(GB.HyperSphere(GB.Point(0.0f0, 0.0f0), 1.0f0), GB.HyperSphere(GB.Point(0.0f0, 0.0f0), 1.0f0), PE2D.Manifold(2.0f0, GB.Vec(1.0f0, 0.0f0))),
                             (GB.HyperSphere(GB.Point(0.0f0, 0.0f0), 1.0f0), GB.HyperSphere(GB.Point(1.0f0, 0.0f0), 1.0f0), PE2D.Manifold(1.0f0, GB.Vec(1.0f0, 0.0f0))),
                             (GB.HyperSphere(GB.Point(0.0f0, 0.0f0), 1.0f0), GB.HyperSphere(GB.Point(2.0f0, 0.0f0), 1.0f0), PE2D.Manifold(0.0f0, GB.Vec(1.0f0, 0.0f0)))]
            test_manifold_list(manifold_list)
        end

        @testset "Rect2D vs. Circle" begin
            manifold_list = [(GB.HyperRectangle(0.0f0, -1.0f0, 5.0f0, 2.0f0), GB.HyperSphere(GB.Point(0.0f0, 0.0f0), 1.0f0), PE2D.Manifold(1.0f0, GB.Vec(-1.0f0, 0.0f0))),
                             (GB.HyperRectangle(0.5f0, -0.5f0, 1.0f0, 1.0f0), GB.HyperSphere(GB.Point(0.0f0, 0.0f0), 1.0f0), PE2D.Manifold(0.5f0, GB.Vec(-1.0f0, 0.0f0))),
                             (GB.HyperRectangle(0.5f0, 0.5f0, 1.0f0, 1.0f0), GB.HyperSphere(GB.Point(0.0f0, 0.0f0), 1.0f0), PE2D.Manifold(1.0f0 - Float32(0.5 * sqrt(2)), GB.Vec(Float32(-0.5 * sqrt(2)), Float32(-0.5 * sqrt(2)))))]
            test_manifold_list(manifold_list)
        end

        @testset "Rect2D vs. Rect2D" begin
            manifold_list = [(GB.HyperRectangle(1.0f0, 2.0f0, 3.0f0, 4.0f0), GB.HyperRectangle(0.0f0, 0.0f0, 1.1f0, 5.0f0), PE2D.Manifold(0.1f0, GB.Vec(-1.0f0, 0.0f0))),
                             (GB.HyperRectangle(1.0f0, 2.0f0, 3.0f0, 4.0f0), GB.HyperRectangle(2.0f0, 0.0f0, 5.0f0, 2.1f0), PE2D.Manifold(0.1f0, GB.Vec(0.0f0, -1.0f0))),
                             (GB.HyperRectangle(1.0f0, 2.0f0, 3.0f0, 4.0f0), GB.HyperRectangle(2.0f0, 0.0f0, 5.0f0, 2.0f0), PE2D.Manifold(0.0f0, GB.Vec(0.0f0, -1.0f0)))]
            test_manifold_list(manifold_list)
        end
    end

    @testset "World simulation" begin
        T = Float32
        NUM_ITER = 500
        DT = 0.01
        FRAME_RATE = 1 / DT

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
        run(world, 500, 0.01)
    end

end
