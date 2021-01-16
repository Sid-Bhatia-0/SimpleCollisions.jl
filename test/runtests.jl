using PhysicsEngine2D
using GeometryBasics
using Test

function test_collision_list(collision_list)
    for (a, b, value) in collision_list
        @test PE2D.is_colliding(a, b) == value
        @test PE2D.is_colliding(b, a) == value
    end
end

@testset "PhysicsEngine2D.jl" begin
    @testset "area computation" begin
        @testset "Rect2D" begin
            a = Rect(1, 2, 3, 4)
            @test area(a) == 12
        end

        @testset "Circle" begin
            a = HyperSphere(Point(0, 0), 1)
            @test area(a) ≈ π
        end
    end

    @testset "collision detection" begin
        @testset "Circle vs. Point" begin
            collision_list = [(HyperSphere(Point(0, 0), 1), Point(0, 0), true),
                              (HyperSphere(Point(0, 0), 1), Point(1, 0), true),
                              (HyperSphere(Point(0, 0), 1), Point(2, 0), false)]
            test_collision_list(collision_list)
        end

        @testset "Circle vs. Circle" begin
            collision_list = [(HyperSphere(Point(0, 0), 1), HyperSphere(Point(0, 0), 1), true),
                              (HyperSphere(Point(0, 0), 1), HyperSphere(Point(1, 0), 1), true),
                              (HyperSphere(Point(0, 0), 1), HyperSphere(Point(2, 0), 1), true),
                              (HyperSphere(Point(0, 0), 1), HyperSphere(Point(3, 0), 1), false)]
            test_collision_list(collision_list)
        end

        @testset "Rect2D vs. Circle" begin
            collision_list = [(Rect(1, 2, 5, 6), HyperSphere(Point(3, 3), 1), true),
                              (Rect(1, 2, 5, 6), HyperSphere(Point(4, 4), 1), true),
                              (Rect(1, 2, 5, 6), HyperSphere(Point(4, 4), 10), true),
                              (Rect(1, 2, 5, 6), HyperSphere(Point(0, 0), 3), true),
                              (Rect(0, 1, 1, 1), HyperSphere(Point(0, 0), 1), true),
                              (Rect(1, 2, 5, 6), HyperSphere(Point(0, 0), 1), false)]
            test_collision_list(collision_list)
        end

        @testset "Rect2D vs. Rect2D" begin
            collision_list = [(Rect(1, 2, 3, 4), Rect(3, 4, 5, 6), true),
                              (Rect(1, 2, 3, 4), Rect(4, 6, 1, 2), true),
                              (Rect(1, 2, 3, 4), Rect(0, 0, 6, 6), true),
                              (Rect(1, 2, 3, 4), Rect(4, 2, 1, 2), true),
                              (Rect(1, 2, 3, 4), Rect(5, 6, 7, 8), false)]
            test_collision_list(collision_list)
        end
    end

    @testset "RigidBody instantiation" begin
        body = PE2D.RigidBody{Float64}()
        body = PE2D.RigidBody{Float32}()
        body = PE2D.RigidBody{Float16}()
    end

    @testset "World instantiation" begin
        world = PE2D.World([])
    end

    @testset "World run" begin
        body = PE2D.RigidBody{Float32}()
        PE2D.set_velocity!(body, Vec2(1.0f0, 0.0f0))

        world = PE2D.World([body])
        run(world, 5, 0.2)
    end
end
