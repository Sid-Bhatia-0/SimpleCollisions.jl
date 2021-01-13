using PhysicsEngine2D
using GeometryBasics
using Test

@testset "PhysicsEngine2D.jl" begin
    @testset "Rect2D intersection" begin
        a = Rect2D(1, 2, 3, 4)
        b = Rect2D(4, 6, 1, 2)
        @test PE2D.is_intersecting(a, b) == true

        a = Rect2D(4, 6, 1, 2)
        b = Rect2D(1, 2, 3, 4)
        @test PE2D.is_intersecting(a, b) == true

        a = Rect2D(1, 2, 3, 4)
        b = Rect2D(4, 2, 1, 2)
        @test PE2D.is_intersecting(a, b) == true

        a = Rect2D(1, 2, 3, 4)
        b = Rect2D(5, 6, 7, 8)
        @test PE2D.is_intersecting(a, b) == false

        a = Rect2D(1, 2, 3, 4)
        b = Rect2D(-2, -2, 1, 2)
        @test PE2D.is_intersecting(a, b) == false
    end

    @testset "Circle intersection" begin
        a = Circle(Point2(0, 0), 1)
        b = Circle(Point2(0, 0), 1)
        @test PE2D.is_intersecting(a, b) == true

        a = Circle(Point2(0, 0), 1)
        b = Circle(Point2(2, 0), 1)
        @test PE2D.is_intersecting(a, b) == true

        a = Circle(Point2(0, 0), 1)
        b = Circle(Point2(3, 0), 1)
        @test PE2D.is_intersecting(a, b) == false
    end
end
