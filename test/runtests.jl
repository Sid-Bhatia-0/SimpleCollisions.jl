using PhysicsEngine2D
using GeometryBasics
using Test

@testset "PhysicsEngine2D.jl" begin
    @testset "Rect2D area" begin
        a = Rect2D(1, 2, 3, 4)
        @test area(a) == 12
    end

    @testset "Circle area" begin
        a = Circle(Point2(0, 0), 1)
        @test area(a) ≈ π
    end

    @testset "Rect2D vs. Rect2D collision" begin
        a = Rect2D(1, 2, 3, 4)
        b = Rect2D(3, 4, 5, 6)
        @test PE2D.is_colliding(a, b) == true

        a = Rect2D(1, 2, 3, 4)
        b = Rect2D(4, 6, 1, 2)
        @test PE2D.is_colliding(a, b) == true

        a = Rect2D(1, 2, 3, 4)
        b = Rect2D(0, 0, 6, 6)
        @test PE2D.is_colliding(a, b) == true

        a = Rect2D(1, 2, 3, 4)
        b = Rect2D(4, 2, 1, 2)
        @test PE2D.is_colliding(a, b) == true

        a = Rect2D(1, 2, 3, 4)
        b = Rect2D(5, 6, 7, 8)
        @test PE2D.is_colliding(a, b) == false
    end

    @testset "Circle vs. Circle collision" begin
        a = Circle(Point2(0, 0), 1)
        b = Circle(Point2(0, 0), 1)
        @test PE2D.is_colliding(a, b) == true

        a = Circle(Point2(0, 0), 1)
        b = Circle(Point2(2, 0), 1)
        @test PE2D.is_colliding(a, b) == true

        a = Circle(Point2(0, 0), 1)
        b = Circle(Point2(3, 0), 1)
        @test PE2D.is_colliding(a, b) == false
    end

    @testset "Circle vs. Point collision" begin
        a = Circle(Point2(0, 0), 1)
        b = Point(0, 0)
        @test PE2D.is_colliding(a, b) == true

        a = Circle(Point2(0, 0), 1)
        b = Point(1, 0)
        @test PE2D.is_colliding(a, b) == true

        a = Circle(Point2(0, 0), 1)
        b = Point(2, 0)
        @test PE2D.is_colliding(a, b) == false
    end

    @testset "Rect2D vs. Circle collision" begin
        a = Rect2D(1, 2, 5, 6)
        b = Circle(Point2(3, 3), 1)
        @test PE2D.is_colliding(a, b) == true

        a = Rect2D(1, 2, 5, 6)
        b = Circle(Point2(4, 4), 1)
        @test PE2D.is_colliding(a, b) == true

        a = Rect2D(1, 2, 5, 6)
        b = Circle(Point2(4, 4), 10)
        @test PE2D.is_colliding(a, b) == true

        a = Rect2D(1, 2, 5, 6)
        b = Circle(Point2(0, 0), 3)
        @test PE2D.is_colliding(a, b) == true

        a = Rect2D(1, 2, 5, 6)
        b = Circle(Point2(0, 0), 1)
        @test PE2D.is_colliding(a, b) == false
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
