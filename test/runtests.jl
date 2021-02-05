import PhysicsEngine2D
import PhysicsEngine2D: PE2D
import GeometryBasics
const GB = GeometryBasics
import LinearAlgebra
const LA = LinearAlgebra
using Test

function test_collision_list(collision_list)
    for (a, b, pos_ba, axes_ba, value) in collision_list
        @test PE2D.is_colliding(a, b, pos_ba, axes_ba) == value
    end
end

function test_manifold_list(manifold_list)
    for (a, b, pos_ba, axes_ba, value) in manifold_list
        manifold_ba = PE2D.Manifold(a, b, pos_ba, axes_ba)

        @show a
        @show b
        @show pos_ba
        @show axes_ba
        @show manifold_ba

        @test PE2D.get_penetration(manifold_ba) ≈ PE2D.get_penetration(value)
        @test PE2D.get_normal(manifold_ba) ≈ PE2D.get_normal(value)
        @test PE2D.get_tangent(manifold_ba) ≈ PE2D.get_tangent(value)
        @test PE2D.get_contact(manifold_ba) ≈ PE2D.get_contact(value)
    end
end

@testset "PhysicsEngine2D.jl" begin
    T = Float32
    penetration = convert(T, 0.01)
    d = convert(T, 0.01)

    origin = zero(GB.Vec2{T})
    VecType = typeof(origin)

    std_axes = PE2D.Axes{T}()
    i_cap = PE2D.get_x_cap(std_axes)
    j_cap = PE2D.get_y_cap(std_axes)

    theta = convert(T, π / 6)
    rotated_axes = PE2D.Axes(theta)

    l1 = GB.Line(convert(GB.Point, -i_cap), convert(GB.Point, i_cap))
    p1_l1 = PE2D.get_point(l1, 1)
    p2_l1 = PE2D.get_point(l1, 2)

    l2 = GB.Line(convert(GB.Point, -2 .* i_cap), convert(GB.Point, 2 .* i_cap))
    p1_l2 = PE2D.get_point(l2, 1)
    p2_l2 = PE2D.get_point(l2, 2)

    c1 = GB.HyperSphere(convert(GB.Point, origin), one(T))
    r_c1 = c1.r

    c2 = GB.HyperSphere(convert(GB.Point, origin), 2 * one(T))
    r_c2 = c2.r

    r1 = GB.Rect(origin .- VecType(1, 0.5), VecType(2, 1))
    top_right_r1 = PE2D.get_top_right(r1)
    half_width_r1 = top_right_r1[1]
    half_height_r1 = top_right_r1[2]
    theta_r1 = atan(half_height_r1, half_width_r1)

    r2 = GB.Rect(origin .- VecType(2, 1), VecType(4, 2))
    top_right_r2 = PE2D.get_top_right(r2)
    half_width_r2 = top_right_r2[1]
    half_height_r2 = top_right_r2[2]
    theta_r2 = atan(half_height_r2, half_width_r2)

    @testset "Area computation" begin
        @testset "Rect2D" begin
            @test PE2D.get_area(GB.Rect(1, 2, 3, 4)) == 12
        end

        @testset "Circle" begin
            @test PE2D.get_area(GB.HyperSphere(GB.Point(0, 0), 1)) ≈ π
        end
    end

    @testset "Collision detection" begin
        @testset "Point2 vs. Point2" begin
            collision_list = [
            # std_axes
            (origin, origin, origin, std_axes, true),
            (origin, origin, i_cap, std_axes, false),
            (origin, origin, j_cap, std_axes, false),

            # rotated_axes
            (origin, origin, origin, rotated_axes, true),
            (origin, origin, i_cap, rotated_axes, false),
            (origin, origin, j_cap, rotated_axes, false),
            ]

            test_collision_list(collision_list)
        end

        @testset "Line segment vs. Point2" begin
            collision_list = [
            # std_axes
            (l1, origin, origin, std_axes, true),

            (l1, origin, -convert(T, 1.01) .* i_cap, std_axes, false),
            (l1, origin, -i_cap, std_axes, true),
            (l1, origin, -convert(T, 0.99) .* i_cap, std_axes, true),
            (l1, origin, convert(T, 0.99) .* i_cap, std_axes, true),
            (l1, origin, i_cap, std_axes, true),
            (l1, origin, convert(T, 1.01) .* i_cap, std_axes, false),

            (l1, origin, -j_cap, std_axes, false),
            (l1, origin, -j_cap ./ 1000, std_axes, false),
            (l1, origin, j_cap, std_axes, false),
            (l1, origin, j_cap ./ 1000, std_axes, false),

            (l1, origin, i_cap .+ j_cap, std_axes, false),

            # rotated_axes
            (l1, origin, origin, rotated_axes, true),

            (l1, origin, -convert(T, 1.01) .* i_cap, rotated_axes, false),
            (l1, origin, -i_cap, rotated_axes, true),
            (l1, origin, -convert(T, 0.99) .* i_cap, rotated_axes, true),
            (l1, origin, convert(T, 0.99) .* i_cap, rotated_axes, true),
            (l1, origin, i_cap, rotated_axes, true),
            (l1, origin, convert(T, 1.01) .* i_cap, rotated_axes, false),

            (l1, origin, -j_cap, rotated_axes, false),
            (l1, origin, -j_cap ./ 1000, rotated_axes, false),
            (l1, origin, j_cap, rotated_axes, false),
            (l1, origin, j_cap ./ 1000, rotated_axes, false),

            (l1, origin, i_cap .+ j_cap, rotated_axes, false),

            # reverse check with std_axes
            (origin, l1, origin, std_axes, true),

            (origin, l1, -convert(T, 1.01) .* i_cap, std_axes, false),
            (origin, l1, -i_cap, std_axes, true),
            (origin, l1, -convert(T, 0.99) .* i_cap, std_axes, true),
            (origin, l1, convert(T, 0.99) .* i_cap, std_axes, true),
            (origin, l1, i_cap, std_axes, true),
            (origin, l1, convert(T, 1.01) .* i_cap, std_axes, false),

            (origin, l1, -j_cap, std_axes, false),
            (origin, l1, -j_cap ./ 1000, std_axes, false),
            (origin, l1, j_cap, std_axes, false),
            (origin, l1, j_cap ./ 1000, std_axes, false),

            (origin, l1, i_cap .+ j_cap, std_axes, false),

            # reverse check with rotated_axes
            (origin, l1, origin, rotated_axes, true),

            (origin, l1, -convert(T, 1.01) .* i_cap, rotated_axes, false),
            (origin, l1, -i_cap, rotated_axes, false),
            (origin, l1, -convert(T, 0.99) .* i_cap, rotated_axes, false),
            (origin, l1, convert(T, 0.99) .* i_cap, rotated_axes, false),
            (origin, l1, i_cap, rotated_axes, false),
            (origin, l1, convert(T, 1.01) .* i_cap, rotated_axes, false),

            (origin, l1, -j_cap, rotated_axes, false),
            (origin, l1, -j_cap ./ 1000, rotated_axes, false),
            (origin, l1, j_cap, rotated_axes, false),
            (origin, l1, j_cap ./ 1000, rotated_axes, false),

            (origin, l1, i_cap .+ j_cap, rotated_axes, false),
            ]

            test_collision_list(collision_list)
        end

        @testset "Line segment vs. Line segment" begin
            collision_list = [
            # std_axes
            (l1, l2, origin, std_axes, true),

            (l1, l2, -convert(T, 2.99) .* i_cap, std_axes, true),
            (l1, l2, -convert(T, 3) .* i_cap, std_axes, true),
            (l1, l2, -convert(T, 3.01) .* i_cap, std_axes, false),
            (l1, l2, convert(T, 2.99) .* i_cap, std_axes, true),
            (l1, l2, convert(T, 3) .* i_cap, std_axes, true),
            (l1, l2, convert(T, 3.01) .* i_cap, std_axes, false),

            (l1, l2, -j_cap, std_axes, false),
            (l1, l2, -j_cap ./ 1000, std_axes, false),
            (l1, l2, j_cap, std_axes, false),
            (l1, l2, j_cap ./ 1000, std_axes, false),

            (l1, l2, i_cap .+ j_cap, std_axes, false),

            # rotated_axes
            (l1, l2, origin, rotated_axes, true),

            (l1, l2, -convert(T, 2.99) .* i_cap, rotated_axes, false),
            (l1, l2, -convert(T, 3) .* i_cap, rotated_axes, false),
            (l1, l2, -convert(T, 3.01) .* i_cap, rotated_axes, false),
            (l1, l2, -convert(T, 1.01) .* i_cap, rotated_axes, false),
            (l1, l2, convert(T, 1.01) .* i_cap, rotated_axes, false),
            (l1, l2, convert(T, 2.99) .* i_cap, rotated_axes, false),
            (l1, l2, convert(T, 3) .* i_cap, rotated_axes, false),
            (l1, l2, convert(T, 3.01) .* i_cap, rotated_axes, false),

            (l1, l2, -j_cap, rotated_axes, false),
            (l1, l2, -j_cap ./ 1000, rotated_axes, true),
            (l1, l2, j_cap, rotated_axes, false),
            (l1, l2, j_cap ./ 1000, rotated_axes, true),

            (l1, l2, i_cap .+ j_cap, rotated_axes, true),
            (l1, l2, i_cap .+ convert(T, 1.01) .* j_cap, rotated_axes, false),
            ]

            test_collision_list(collision_list)
        end

        @testset "Circle vs. Point2" begin
            collision_list = [
            # std_axes
            (c1, origin, origin, std_axes, true),

            (c1, origin, -convert(T, 1.01) .* i_cap, std_axes, false),
            (c1, origin, -i_cap, std_axes, true),
            (c1, origin, -convert(T, 0.99) .* i_cap, std_axes, true),
            (c1, origin, convert(T, 0.99) .* i_cap, std_axes, true),
            (c1, origin, i_cap, std_axes, true),
            (c1, origin, convert(T, 1.01) .* i_cap, std_axes, false),

            (c1, origin, -convert(T, 1.01) .* j_cap, std_axes, false),
            (c1, origin, -j_cap, std_axes, true),
            (c1, origin, -convert(T, 0.99) .* j_cap, std_axes, true),
            (c1, origin, convert(T, 0.99) .* j_cap, std_axes, true),
            (c1, origin, j_cap, std_axes, true),
            (c1, origin, convert(T, 1.01) .* j_cap, std_axes, false),

            (c1, origin, i_cap .+ j_cap, std_axes, false),
            (c1, origin, (i_cap .+ j_cap) ./ 2, std_axes, true),

            # rotated_axes
            (c1, origin, origin, rotated_axes, true),

            (c1, origin, -convert(T, 1.01) .* i_cap, rotated_axes, false),
            (c1, origin, -i_cap, rotated_axes, true),
            (c1, origin, -convert(T, 0.99) .* i_cap, rotated_axes, true),
            (c1, origin, convert(T, 0.99) .* i_cap, rotated_axes, true),
            (c1, origin, i_cap, rotated_axes, true),
            (c1, origin, convert(T, 1.01) .* i_cap, rotated_axes, false),

            (c1, origin, -convert(T, 1.01) .* j_cap, rotated_axes, false),
            (c1, origin, -j_cap, rotated_axes, true),
            (c1, origin, -convert(T, 0.99) .* j_cap, rotated_axes, true),
            (c1, origin, convert(T, 0.99) .* j_cap, rotated_axes, true),
            (c1, origin, j_cap, rotated_axes, true),
            (c1, origin, convert(T, 1.01) .* j_cap, rotated_axes, false),

            (c1, origin, i_cap .+ j_cap, rotated_axes, false),
            (c1, origin, (i_cap .+ j_cap) ./ 2, rotated_axes, true),

            # reverse check with std_axes
            (origin, c1, origin, std_axes, true),

            (origin, c1, -convert(T, 1.01) .* i_cap, std_axes, false),
            (origin, c1, -i_cap, std_axes, true),
            (origin, c1, -convert(T, 0.99) .* i_cap, std_axes, true),
            (origin, c1, convert(T, 0.99) .* i_cap, std_axes, true),
            (origin, c1, i_cap, std_axes, true),
            (origin, c1, convert(T, 1.01) .* i_cap, std_axes, false),

            (origin, c1, -convert(T, 1.01) .* j_cap, std_axes, false),
            (origin, c1, -j_cap, std_axes, true),
            (origin, c1, -convert(T, 0.99) .* j_cap, std_axes, true),
            (origin, c1, convert(T, 0.99) .* j_cap, std_axes, true),
            (origin, c1, j_cap, std_axes, true),
            (origin, c1, convert(T, 1.01) .* j_cap, std_axes, false),

            (origin, c1, i_cap .+ j_cap, std_axes, false),
            (origin, c1, (i_cap .+ j_cap) ./ 2, std_axes, true),

            # reverse check with rotated_axes
            (origin, c1, origin, rotated_axes, true),

            (origin, c1, -convert(T, 1.01) .* i_cap, rotated_axes, false),
            (origin, c1, -i_cap, rotated_axes, true),
            (origin, c1, -convert(T, 0.99) .* i_cap, rotated_axes, true),
            (origin, c1, convert(T, 0.99) .* i_cap, rotated_axes, true),
            (origin, c1, i_cap, rotated_axes, true),
            (origin, c1, convert(T, 1.01) .* i_cap, rotated_axes, false),

            (origin, c1, -convert(T, 1.01) .* j_cap, rotated_axes, false),
            (origin, c1, -j_cap, rotated_axes, true),
            (origin, c1, -convert(T, 0.99) .* j_cap, rotated_axes, true),
            (origin, c1, convert(T, 0.99) .* j_cap, rotated_axes, true),
            (origin, c1, j_cap, rotated_axes, true),
            (origin, c1, convert(T, 1.01) .* j_cap, rotated_axes, false),

            (origin, c1, i_cap .+ j_cap, rotated_axes, false),
            (origin, c1, (i_cap .+ j_cap) ./ 2, rotated_axes, true),
            ]

            test_collision_list(collision_list)
        end

        @testset "Circle vs. Line segment" begin
            collision_list = [
            # std_axes
            (l1, c2, origin, std_axes, true),
            (l1, c2, origin, std_axes, true),

            (l1, c2, -convert(T, 3.01) .* i_cap, std_axes, false),
            (l1, c2, -convert(T, 3) .* i_cap, std_axes, true),
            (l1, c2, -convert(T, 2.99) .* i_cap, std_axes, true),
            (l1, c2, -convert(T, 1.01) .* i_cap, std_axes, true),
            (l1, c2, convert(T, 1.01) .* i_cap, std_axes, true),
            (l1, c2, convert(T, 2.99) .* i_cap, std_axes, true),
            (l1, c2, convert(T, 3) .* i_cap, std_axes, true),
            (l1, c2, convert(T, 3.01) .* i_cap, std_axes, false),

            (l1, c2, -convert(T, 2.01) .* j_cap, std_axes, false),
            (l1, c2, -convert(T, 2) .* j_cap, std_axes, true),
            (l1, c2, -convert(T, 1.99) .* j_cap, std_axes, true),
            (l1, c2, -convert(T, 1.01) .* j_cap, std_axes, true),
            (l1, c2, convert(T, 1.01) .* j_cap, std_axes, true),
            (l1, c2, convert(T, 1.99) .* j_cap, std_axes, true),
            (l1, c2, convert(T, 2) .* j_cap, std_axes, true),
            (l1, c2, convert(T, 2.01) .* j_cap, std_axes, false),

            (l1, c2, i_cap .* j_cap, std_axes, true),

            # rotated_axes
            (l1, c2, origin, rotated_axes, true),
            (l1, c2, origin, rotated_axes, true),

            (l1, c2, -convert(T, 3.01) .* i_cap, rotated_axes, false),
            (l1, c2, -convert(T, 3) .* i_cap, rotated_axes, true),
            (l1, c2, -convert(T, 2.99) .* i_cap, rotated_axes, true),
            (l1, c2, -convert(T, 1.01) .* i_cap, rotated_axes, true),
            (l1, c2, convert(T, 1.01) .* i_cap, rotated_axes, true),
            (l1, c2, convert(T, 2.99) .* i_cap, rotated_axes, true),
            (l1, c2, convert(T, 3) .* i_cap, rotated_axes, true),
            (l1, c2, convert(T, 3.01) .* i_cap, rotated_axes, false),

            (l1, c2, -convert(T, 2.01) .* j_cap, rotated_axes, false),
            (l1, c2, -convert(T, 2) .* j_cap, rotated_axes, true),
            (l1, c2, -convert(T, 1.99) .* j_cap, rotated_axes, true),
            (l1, c2, -convert(T, 1.01) .* j_cap, rotated_axes, true),
            (l1, c2, convert(T, 1.01) .* j_cap, rotated_axes, true),
            (l1, c2, convert(T, 1.99) .* j_cap, rotated_axes, true),
            (l1, c2, convert(T, 2) .* j_cap, rotated_axes, true),
            (l1, c2, convert(T, 2.01) .* j_cap, rotated_axes, false),

            (l1, c2, i_cap .* j_cap, rotated_axes, true),

            # reverse check with std_axes
            (c2, l1, origin, std_axes, true),
            (c2, l1, origin, std_axes, true),

            (c2, l1, -convert(T, 3.01) .* i_cap, std_axes, false),
            (c2, l1, -convert(T, 3) .* i_cap, std_axes, true),
            (c2, l1, -convert(T, 2.99) .* i_cap, std_axes, true),
            (c2, l1, -convert(T, 1.01) .* i_cap, std_axes, true),
            (c2, l1, convert(T, 1.01) .* i_cap, std_axes, true),
            (c2, l1, convert(T, 2.99) .* i_cap, std_axes, true),
            (c2, l1, convert(T, 3) .* i_cap, std_axes, true),
            (c2, l1, convert(T, 3.01) .* i_cap, std_axes, false),

            (c2, l1, -convert(T, 2.01) .* j_cap, std_axes, false),
            (c2, l1, -convert(T, 2) .* j_cap, std_axes, true),
            (c2, l1, -convert(T, 1.99) .* j_cap, std_axes, true),
            (c2, l1, -convert(T, 1.01) .* j_cap, std_axes, true),
            (c2, l1, convert(T, 1.01) .* j_cap, std_axes, true),
            (c2, l1, convert(T, 1.99) .* j_cap, std_axes, true),
            (c2, l1, convert(T, 2) .* j_cap, std_axes, true),
            (c2, l1, convert(T, 2.01) .* j_cap, std_axes, false),

            (c2, l1, i_cap .* j_cap, std_axes, true),

            # reverse check with rotated_axes
            (c2, l1, origin, rotated_axes, true),
            (c2, l1, origin, rotated_axes, true),

            (c2, l1, -convert(T, 3.01) .* i_cap, rotated_axes, false),
            (c2, l1, -convert(T, 3) .* i_cap, rotated_axes, false),
            (c2, l1, convert(T, 3) .* i_cap, rotated_axes, false),
            (c2, l1, convert(T, 3.01) .* i_cap, rotated_axes, false),

            (c2, l1, -convert(T, 2.01) .* j_cap, rotated_axes, true),
            (c2, l1, -convert(T, 2) .* j_cap, rotated_axes, true),
            (c2, l1, -convert(T, 1.99) .* j_cap, rotated_axes, true),
            (c2, l1, -convert(T, 1.01) .* j_cap, rotated_axes, true),
            (c2, l1, convert(T, 1.01) .* j_cap, rotated_axes, true),
            (c2, l1, convert(T, 1.99) .* j_cap, rotated_axes, true),
            (c2, l1, convert(T, 2) .* j_cap, rotated_axes, true),
            (c2, l1, convert(T, 2.01) .* j_cap, rotated_axes, true),

            (c2, l1, convert(T, 1.5 * sqrt(3) - 0.01) .* i_cap .+ convert(T, 1.5 - 0.01) .* j_cap, rotated_axes, true),
            (c2, l1, convert(T, 1.5 * sqrt(3)) .* i_cap .+ convert(T, 1.5) .* j_cap, rotated_axes, true),
            (c2, l1, convert(T, 1.5 * sqrt(3) + 0.01) .* i_cap .+ convert(T, 1.5 + 0.01) .* j_cap, rotated_axes, false),
            ]

            test_collision_list(collision_list)
        end

        @testset "Circle vs. Circle" begin
            collision_list = [
            # std_axes
            (c1, c2, origin, std_axes, true),
            (c1, c2, convert(T, 2) .* i_cap, std_axes, true),

            (c1, c2, -convert(T, 3.01) .* i_cap, std_axes, false),
            (c1, c2, -convert(T, 3) .* i_cap, std_axes, true),
            (c1, c2, -convert(T, 2.99) .* i_cap, std_axes, true),
            (c1, c2, -convert(T, 1.01) .* i_cap, std_axes, true),
            (c1, c2, convert(T, 1.01) .* i_cap, std_axes, true),
            (c1, c2, convert(T, 2.99) .* i_cap, std_axes, true),
            (c1, c2, convert(T, 3) .* i_cap, std_axes, true),
            (c1, c2, convert(T, 3.01) .* i_cap, std_axes, false),

            (c1, c2, -convert(T, 3.01) .* j_cap, std_axes, false),
            (c1, c2, -convert(T, 3) .* j_cap, std_axes, true),
            (c1, c2, -convert(T, 2.99) .* j_cap, std_axes, true),
            (c1, c2, -convert(T, 1.01) .* j_cap, std_axes, true),
            (c1, c2, convert(T, 1.01) .* j_cap, std_axes, true),
            (c1, c2, convert(T, 2.99) .* j_cap, std_axes, true),
            (c1, c2, convert(T, 3) .* j_cap, std_axes, true),
            (c1, c2, convert(T, 3.01) .* j_cap, std_axes, false),

            (c1, c2, convert(T, 3 / sqrt(2) - 0.01) .* i_cap .+ convert(T, 3 / sqrt(2) - 0.01) .* j_cap, std_axes, true),
            (c1, c2, convert(T, 3 / sqrt(2)) .* i_cap .+ convert(T, 3 / sqrt(2)) .* j_cap, std_axes, true),
            (c1, c2, convert(T, 3 / sqrt(2) + 0.01) .* i_cap .+ convert(T, 3 / sqrt(2) + 0.01) .* j_cap, std_axes, false),

            # rotated_axes
            (c1, c2, origin, rotated_axes, true),
            (c1, c2, convert(T, 2) .* i_cap, rotated_axes, true),

            (c1, c2, -convert(T, 3.01) .* i_cap, rotated_axes, false),
            (c1, c2, -convert(T, 3) .* i_cap, rotated_axes, true),
            (c1, c2, -convert(T, 2.99) .* i_cap, rotated_axes, true),
            (c1, c2, -convert(T, 1.01) .* i_cap, rotated_axes, true),
            (c1, c2, convert(T, 1.01) .* i_cap, rotated_axes, true),
            (c1, c2, convert(T, 2.99) .* i_cap, rotated_axes, true),
            (c1, c2, convert(T, 3) .* i_cap, rotated_axes, true),
            (c1, c2, convert(T, 3.01) .* i_cap, rotated_axes, false),

            (c1, c2, -convert(T, 3.01) .* j_cap, rotated_axes, false),
            (c1, c2, -convert(T, 3) .* j_cap, rotated_axes, true),
            (c1, c2, -convert(T, 2.99) .* j_cap, rotated_axes, true),
            (c1, c2, -convert(T, 1.01) .* j_cap, rotated_axes, true),
            (c1, c2, convert(T, 1.01) .* j_cap, rotated_axes, true),
            (c1, c2, convert(T, 2.99) .* j_cap, rotated_axes, true),
            (c1, c2, convert(T, 3) .* j_cap, rotated_axes, true),
            (c1, c2, convert(T, 3.01) .* j_cap, rotated_axes, false),

            (c1, c2, convert(T, 3 / sqrt(2) - 0.01) .* i_cap .+ convert(T, 3 / sqrt(2) - 0.01) .* j_cap, rotated_axes, true),
            (c1, c2, convert(T, 3 / sqrt(2)) .* i_cap .+ convert(T, 3 / sqrt(2)) .* j_cap, rotated_axes, true),
            (c1, c2, convert(T, 3 / sqrt(2) + 0.01) .* i_cap .+ convert(T, 3 / sqrt(2) + 0.01) .* j_cap, rotated_axes, false),
            ]

            test_collision_list(collision_list)
        end

        @testset "Rect2D vs. Point2" begin
            collision_list = [
            # std_axes
            (r1, origin, origin, std_axes, true),

            (r1, origin, -convert(T, 1.01) .* i_cap, std_axes, false),
            (r1, origin, -i_cap, std_axes, true),
            (r1, origin, -convert(T, 0.99) .* i_cap, std_axes, true),
            (r1, origin, convert(T, 0.99) .* i_cap, std_axes, true),
            (r1, origin, i_cap, std_axes, true),
            (r1, origin, convert(T, 1.01) .* i_cap, std_axes, false),

            (r1, origin, -convert(T, 0.51) .* j_cap, std_axes, false),
            (r1, origin, -convert(T, 0.50) .* j_cap, std_axes, true),
            (r1, origin, -convert(T, 0.49) .* j_cap, std_axes, true),
            (r1, origin, convert(T, 0.49) .* j_cap, std_axes, true),
            (r1, origin, convert(T, 0.50) .* j_cap, std_axes, true),
            (r1, origin, convert(T, 0.51) .* j_cap, std_axes, false),

            (r1, origin, maximum(r1) .- 0.01, std_axes, true),
            (r1, origin, maximum(r1), std_axes, true),
            (r1, origin, maximum(r1) .+ 0.01, std_axes, false),

            # rotated_axes
            (r1, origin, origin, rotated_axes, true),

            (r1, origin, -convert(T, 1.01) .* i_cap, rotated_axes, false),
            (r1, origin, -i_cap, rotated_axes, true),
            (r1, origin, -convert(T, 0.99) .* i_cap, rotated_axes, true),
            (r1, origin, convert(T, 0.99) .* i_cap, rotated_axes, true),
            (r1, origin, i_cap, rotated_axes, true),
            (r1, origin, convert(T, 1.01) .* i_cap, rotated_axes, false),

            (r1, origin, -convert(T, 0.51) .* j_cap, rotated_axes, false),
            (r1, origin, -convert(T, 0.50) .* j_cap, rotated_axes, true),
            (r1, origin, -convert(T, 0.49) .* j_cap, rotated_axes, true),
            (r1, origin, convert(T, 0.49) .* j_cap, rotated_axes, true),
            (r1, origin, convert(T, 0.50) .* j_cap, rotated_axes, true),
            (r1, origin, convert(T, 0.51) .* j_cap, rotated_axes, false),

            (r1, origin, maximum(r1) .- 0.01, rotated_axes, true),
            (r1, origin, maximum(r1), rotated_axes, true),
            (r1, origin, maximum(r1) .+ 0.01, rotated_axes, false),

            # reverse check with std_axes
            (origin, r1, origin, std_axes, true),

            (origin, r1, -convert(T, 1.01) .* i_cap, std_axes, false),
            (origin, r1, -i_cap, std_axes, true),
            (origin, r1, -convert(T, 0.99) .* i_cap, std_axes, true),
            (origin, r1, convert(T, 0.99) .* i_cap, std_axes, true),
            (origin, r1, i_cap, std_axes, true),
            (origin, r1, convert(T, 1.01) .* i_cap, std_axes, false),

            (origin, r1, -convert(T, 0.51) .* j_cap, std_axes, false),
            (origin, r1, -convert(T, 0.50) .* j_cap, std_axes, true),
            (origin, r1, -convert(T, 0.49) .* j_cap, std_axes, true),
            (origin, r1, convert(T, 0.49) .* j_cap, std_axes, true),
            (origin, r1, convert(T, 0.50) .* j_cap, std_axes, true),
            (origin, r1, convert(T, 0.51) .* j_cap, std_axes, false),

            (origin, r1, maximum(r1) .- 0.01, std_axes, true),
            (origin, r1, maximum(r1), std_axes, true),
            (origin, r1, maximum(r1) .+ 0.01, std_axes, false),

            # reverse check with rotated_axes
            (origin, r1, origin, rotated_axes, true),

            (origin, r1, -convert(T, 1.01) .* i_cap, rotated_axes, false),
            (origin, r1, -i_cap, rotated_axes, true),
            (origin, r1, -convert(T, 0.99) .* i_cap, rotated_axes, true),
            (origin, r1, convert(T, 0.99) .* i_cap, rotated_axes, true),
            (origin, r1, i_cap, rotated_axes, true),
            (origin, r1, convert(T, 1.01) .* i_cap, rotated_axes, false),

            (origin, r1, -convert(T, 1 / sqrt(3) + 0.01) .* j_cap, rotated_axes, false),
            (origin, r1, -convert(T, 1 / sqrt(3)) .* j_cap, rotated_axes, true),
            (origin, r1, -convert(T, 1 / sqrt(3) - 0.01) .* j_cap, rotated_axes, true),
            (origin, r1, convert(T, 1 / sqrt(3) - 0.01) .* j_cap, rotated_axes, true),
            (origin, r1, convert(T, 1 / sqrt(3)) .* j_cap, rotated_axes, true),
            (origin, r1, convert(T, 1 / sqrt(3) + 0.01) .* j_cap, rotated_axes, false),

            (origin, r1, PE2D.rotate(PE2D.get_vertices(r1)[1], rotated_axes) .- 0.01, rotated_axes, false),
            (origin, r1, PE2D.rotate(PE2D.get_vertices(r1)[1], rotated_axes), rotated_axes, true),
            (origin, r1, PE2D.rotate(PE2D.get_vertices(r1)[1], rotated_axes) .+ 0.01, rotated_axes, true),
            ]

            test_collision_list(collision_list)
        end

        @testset "Rect2D vs. Line segment" begin
            collision_list = [
            # std_axes
            (r1, l1, origin, std_axes, true),
            (r1, l2, origin, std_axes, true),

            (r1, l1, -convert(T, 2.01) .* i_cap, std_axes, false),
            (r1, l1, -convert(T, 2) .* i_cap, std_axes, true),
            (r1, l1, -convert(T, 1.99) .* i_cap, std_axes, true),
            (r1, l1, convert(T, 1.99) .* i_cap, std_axes, true),
            (r1, l1, convert(T, 2) .* i_cap, std_axes, true),
            (r1, l1, convert(T, 2.01) .* i_cap, std_axes, false),

            (r1, l1, -convert(T, 0.51) .* j_cap, std_axes, false),
            (r1, l1, -convert(T, 0.50) .* j_cap, std_axes, true),
            (r1, l1, -convert(T, 0.49) .* j_cap, std_axes, true),
            (r1, l1, convert(T, 0.49) .* j_cap, std_axes, true),
            (r1, l1, convert(T, 0.50) .* j_cap, std_axes, true),
            (r1, l1, convert(T, 0.51) .* j_cap, std_axes, false),

            (r1, l1, maximum(r1) .- 0.01, std_axes, true),
            (r1, l1, maximum(r1), std_axes, true),
            (r1, l1, maximum(r1) .+ 0.01, std_axes, false),

            # rotated_axes
            (r1, l1, origin, rotated_axes, true),
            (r1, l2, origin, rotated_axes, true),

            (r1, l1, -convert(T, 2.01) .* i_cap, rotated_axes, false),
            (r1, l1, -convert(T, 2) .* i_cap, rotated_axes, false),
            (r1, l1, -convert(T, 1.99) .* i_cap, rotated_axes, false),
            (r1, l1, convert(T, 2.01) .* i_cap, rotated_axes, false),
            (r1, l1, convert(T, 2) .* i_cap, rotated_axes, false),
            (r1, l1, convert(T, 1.99) .* i_cap, rotated_axes, false),

            (r1, l1, -convert(T, 0.51) .* j_cap, rotated_axes, true),
            (r1, l1, -convert(T, 0.50) .* j_cap, rotated_axes, true),
            (r1, l1, -convert(T, 0.49) .* j_cap, rotated_axes, true),
            (r1, l1, convert(T, 0.49) .* j_cap, rotated_axes, true),
            (r1, l1, convert(T, 0.50) .* j_cap, rotated_axes, true),
            (r1, l1, convert(T, 0.51) .* j_cap, rotated_axes, true),

            (r1, l1, maximum(r1) .- 0.01, rotated_axes, true),
            (r1, l1, maximum(r1), rotated_axes, true),
            (r1, l1, maximum(r1) .+ 0.01, rotated_axes, true),

            # reverse check with std_axes
            (r1, l1, origin, std_axes, true),
            (r1, l2, origin, std_axes, true),

            (r1, l1, -convert(T, 2.01) .* i_cap, std_axes, false),
            (r1, l1, -convert(T, 2) .* i_cap, std_axes, true),
            (r1, l1, -convert(T, 1.99) .* i_cap, std_axes, true),
            (r1, l1, convert(T, 1.99) .* i_cap, std_axes, true),
            (r1, l1, convert(T, 2) .* i_cap, std_axes, true),
            (r1, l1, convert(T, 2.01) .* i_cap, std_axes, false),

            (r1, l1, -convert(T, 0.51) .* j_cap, std_axes, false),
            (r1, l1, -convert(T, 0.50) .* j_cap, std_axes, true),
            (r1, l1, -convert(T, 0.49) .* j_cap, std_axes, true),
            (r1, l1, convert(T, 0.49) .* j_cap, std_axes, true),
            (r1, l1, convert(T, 0.50) .* j_cap, std_axes, true),
            (r1, l1, convert(T, 0.51) .* j_cap, std_axes, false),

            (r1, l1, maximum(r1) .- 0.01, std_axes, true),
            (r1, l1, maximum(r1), std_axes, true),
            (r1, l1, maximum(r1) .+ 0.01, std_axes, false),

            # reverse check with rotated_axes
            (l1, r1, origin, rotated_axes, true),
            (l2, r1, origin, rotated_axes, true),

            (l1, r1, -convert(T, 2.01) .* i_cap, rotated_axes, false),
            (l1, r1, -convert(T, 2) .* i_cap, rotated_axes, true),
            (l1, r1, -convert(T, 1.99) .* i_cap, rotated_axes, true),
            (l1, r1, convert(T, 1.99) .* i_cap, rotated_axes, true),
            (l1, r1, convert(T, 2) .* i_cap, rotated_axes, true),
            (l1, r1, convert(T, 2.01) .* i_cap, rotated_axes, false),

            (l1, r1, -convert(T, 0.51) .* j_cap, rotated_axes, true),
            (l1, r1, -convert(T, 0.50) .* j_cap, rotated_axes, true),
            (l1, r1, -convert(T, 0.49) .* j_cap, rotated_axes, true),
            (l1, r1, convert(T, 0.49) .* j_cap, rotated_axes, true),
            (l1, r1, convert(T, 0.50) .* j_cap, rotated_axes, true),
            (l1, r1, convert(T, 0.51) .* j_cap, rotated_axes, true),

            (l1, r1, maximum(r1) .- 0.01, rotated_axes, true),
            (l1, r1, maximum(r1), rotated_axes, true),
            (l1, r1, maximum(r1) .+ 0.01, rotated_axes, true),
            ]

            test_collision_list(collision_list)
        end

        @testset "Rect2D vs. Circle" begin
            top_right = PE2D.get_top_right(r1)
            theta_r1 = atan(top_right[2], top_right[1])

            collision_list = [
            # std_axes
            (r1, c1, origin, std_axes, true),

            (r1, c1, -convert(T, 2.01) .* i_cap, std_axes, false),
            (r1, c1, -convert(T, 2) .* i_cap, std_axes, true),
            (r1, c1, -convert(T, 1.99) .* i_cap, std_axes, true),
            (r1, c1, convert(T, 1.99) .* i_cap, std_axes, true),
            (r1, c1, convert(T, 2) .* i_cap, std_axes, true),
            (r1, c1, convert(T, 2.01) .* i_cap, std_axes, false),

            (r1, c1, -convert(T, 1.51) .* j_cap, std_axes, false),
            (r1, c1, -convert(T, 1.50) .* j_cap, std_axes, true),
            (r1, c1, -convert(T, 1.49) .* j_cap, std_axes, true),
            (r1, c1, convert(T, 1.49) .* j_cap, std_axes, true),
            (r1, c1, convert(T, 1.50) .* j_cap, std_axes, true),
            (r1, c1, convert(T, 1.51) .* j_cap, std_axes, false),

            (r1, c1, top_right .+ VecType(c1.r * cos(theta_r1), c1.r * sin(theta_r1)) .- 0.01, std_axes, true),
            (r1, c1, top_right .+ VecType(c1.r * cos(theta_r1), c1.r * sin(theta_r1)), std_axes, true),
            (r1, c1, top_right .+ VecType(c1.r * cos(theta_r1), c1.r * sin(theta_r1)) .+ 0.01, std_axes, false),

            # rotated_axes
            (r1, c1, origin, rotated_axes, true),

            (r1, c1, -convert(T, 2.01) .* i_cap, rotated_axes, false),
            (r1, c1, -convert(T, 2) .* i_cap, rotated_axes, true),
            (r1, c1, -convert(T, 1.99) .* i_cap, rotated_axes, true),
            (r1, c1, convert(T, 1.99) .* i_cap, rotated_axes, true),
            (r1, c1, convert(T, 2) .* i_cap, rotated_axes, true),
            (r1, c1, convert(T, 2.01) .* i_cap, rotated_axes, false),

            (r1, c1, -convert(T, 1.51) .* j_cap, rotated_axes, false),
            (r1, c1, -convert(T, 1.50) .* j_cap, rotated_axes, true),
            (r1, c1, -convert(T, 1.49) .* j_cap, rotated_axes, true),
            (r1, c1, convert(T, 1.49) .* j_cap, rotated_axes, true),
            (r1, c1, convert(T, 1.50) .* j_cap, rotated_axes, true),
            (r1, c1, convert(T, 1.51) .* j_cap, rotated_axes, false),

            (r1, c1, top_right .+ VecType(c1.r * cos(theta_r1), c1.r * sin(theta_r1)) .- 0.01, rotated_axes, true),
            (r1, c1, top_right .+ VecType(c1.r * cos(theta_r1), c1.r * sin(theta_r1)), rotated_axes, true),
            (r1, c1, top_right .+ VecType(c1.r * cos(theta_r1), c1.r * sin(theta_r1)) .+ 0.01, rotated_axes, false),

            # reverse check with std_axes
            (c1, r1, origin, std_axes, true),

            (c1, r1, -convert(T, 2.01) .* i_cap, std_axes, false),
            (c1, r1, -convert(T, 2) .* i_cap, std_axes, true),
            (c1, r1, -convert(T, 1.99) .* i_cap, std_axes, true),
            (c1, r1, convert(T, 1.99) .* i_cap, std_axes, true),
            (c1, r1, convert(T, 2) .* i_cap, std_axes, true),
            (c1, r1, convert(T, 2.01) .* i_cap, std_axes, false),

            (c1, r1, -convert(T, 1.51) .* j_cap, std_axes, false),
            (c1, r1, -convert(T, 1.50) .* j_cap, std_axes, true),
            (c1, r1, -convert(T, 1.49) .* j_cap, std_axes, true),
            (c1, r1, convert(T, 1.49) .* j_cap, std_axes, true),
            (c1, r1, convert(T, 1.50) .* j_cap, std_axes, true),
            (c1, r1, convert(T, 1.51) .* j_cap, std_axes, false),

            (c1, r1, top_right .+ VecType(c1.r * cos(theta_r1), c1.r * sin(theta_r1)) .- 0.01, std_axes, true),
            (c1, r1, top_right .+ VecType(c1.r * cos(theta_r1), c1.r * sin(theta_r1)), std_axes, true),
            (c1, r1, top_right .+ VecType(c1.r * cos(theta_r1), c1.r * sin(theta_r1)) .+ 0.01, std_axes, false),

            # reverse check with rotated_axes
            (c1, r1, origin, rotated_axes, true),

            (c1, r1, -convert(T, 2.01) .* i_cap, rotated_axes, true),
            (c1, r1, -convert(T, 2) .* i_cap, rotated_axes, true),
            (c1, r1, -convert(T, 1.99) .* i_cap, rotated_axes, true),
            (c1, r1, convert(T, 1.99) .* i_cap, rotated_axes, true),
            (c1, r1, convert(T, 2) .* i_cap, rotated_axes, true),
            (c1, r1, convert(T, 2.01) .* i_cap, rotated_axes, true),

            (c1, r1, -convert(T, 1.51) .* j_cap, rotated_axes, true),
            (c1, r1, -convert(T, 1.50) .* j_cap, rotated_axes, true),
            (c1, r1, -convert(T, 1.49) .* j_cap, rotated_axes, true),
            (c1, r1, convert(T, 1.49) .* j_cap, rotated_axes, true),
            (c1, r1, convert(T, 1.50) .* j_cap, rotated_axes, true),
            (c1, r1, convert(T, 1.51) .* j_cap, rotated_axes, true),

            (c1, r1, convert(T, 1.99) .* PE2D.get_x_cap(rotated_axes), rotated_axes, true),
            (c1, r1, convert(T, 2) .* PE2D.get_x_cap(rotated_axes), rotated_axes, true),
            (c1, r1, convert(T, 2.01) .* PE2D.get_x_cap(rotated_axes), rotated_axes, false),
            ]

            test_collision_list(collision_list)
        end

        @testset "Rect2D vs. Rect2D" begin
            collision_list = [
            # std_axes
            (r1, r2, origin, std_axes, true),

            (r1, r2, -convert(T, 3.01) .* i_cap, std_axes, false),
            (r1, r2, -convert(T, 3) .* i_cap, std_axes, true),
            (r1, r2, -convert(T, 2.99) .* i_cap, std_axes, true),
            (r1, r2, convert(T, 2.99) .* i_cap, std_axes, true),
            (r1, r2, convert(T, 3) .* i_cap, std_axes, true),
            (r1, r2, convert(T, 3.01) .* i_cap, std_axes, false),

            (r1, r2, -convert(T, 1.51) .* j_cap, std_axes, false),
            (r1, r2, -convert(T, 1.50) .* j_cap, std_axes, true),
            (r1, r2, -convert(T, 1.49) .* j_cap, std_axes, true),
            (r1, r2, convert(T, 1.49) .* j_cap, std_axes, true),
            (r1, r2, convert(T, 1.50) .* j_cap, std_axes, true),
            (r1, r2, convert(T, 1.51) .* j_cap, std_axes, false),

            (r1, r2, maximum(r1) .+ maximum(r2) .- convert(T, 0.01), std_axes, true),
            (r1, r2, maximum(r1) .+ maximum(r2), std_axes, true),
            (r1, r2, maximum(r1) .+ maximum(r2) .+ convert(T, 0.01), std_axes, false),

            # rotated_axes
            (r1, r2, origin, rotated_axes, true),

            (r1, r2, -convert(T, 3.01) .* i_cap, rotated_axes, true),
            (r1, r2, -convert(T, 3) .* i_cap, rotated_axes, true),
            (r1, r2, -convert(T, 2.99) .* i_cap, rotated_axes, true),
            (r1, r2, convert(T, 2.99) .* i_cap, rotated_axes, true),
            (r1, r2, convert(T, 3) .* i_cap, rotated_axes, true),
            (r1, r2, convert(T, 3.01) .* i_cap, rotated_axes, true),

            (r1, r2, -convert(T, 1.51) .* j_cap, rotated_axes, true),
            (r1, r2, -convert(T, 1.50) .* j_cap, rotated_axes, true),
            (r1, r2, -convert(T, 1.49) .* j_cap, rotated_axes, true),
            (r1, r2, convert(T, 1.49) .* j_cap, rotated_axes, true),
            (r1, r2, convert(T, 1.50) .* j_cap, rotated_axes, true),
            (r1, r2, convert(T, 1.51) .* j_cap, rotated_axes, true),

            (r1, r2, VecType(0, half_height_r1 - penetration) .+ LA.norm(top_right_r2) .* VecType(cos(theta + theta_r2), sin(theta + theta_r2)), rotated_axes, true),
            (r1, r2, VecType(0, half_height_r1 + penetration) .+ LA.norm(top_right_r2) .* VecType(cos(theta + theta_r2), sin(theta + theta_r2)), rotated_axes, false),
            ]

            test_collision_list(collision_list)
        end
    end

    @testset "Manifold generation" begin
        @testset "Circle vs. Circle" begin
            manifold_list = [
            # std_axes
            (c1, c2, origin, std_axes, PE2D.Manifold(r_c1 + r_c2, PE2D.rotate_minus_90(std_axes), (r_c1 - (r_c1 + r_c2) / 2) .* i_cap)),

            (c1, c2, r_c1 * i_cap, std_axes, PE2D.Manifold(r_c2, PE2D.rotate_minus_90(std_axes), (r_c1 - r_c2 / 2) .* i_cap)),
            (c1, c2, r_c2 * i_cap, std_axes, PE2D.Manifold(r_c1, PE2D.rotate_minus_90(std_axes), (r_c1 - r_c1 / 2) .* i_cap)),

            (c1, c2, r_c1 * j_cap, std_axes, PE2D.Manifold(r_c2, std_axes, (r_c1 - r_c2 / 2) .* j_cap)),
            (c1, c2, r_c2 * j_cap, std_axes, PE2D.Manifold(r_c1, std_axes, (r_c1 - r_c1 / 2) .* j_cap)),

            (c1, c2, r_c1 .* (i_cap .+ j_cap) ./ convert(T, sqrt(2)), std_axes, PE2D.Manifold(r_c2, PE2D.Axes(convert(T, -pi / 4)), (r_c1 - r_c2 / 2) .* (i_cap .+ j_cap) ./ convert(T, sqrt(2)))),
            (c1, c2, r_c2 .* (i_cap .+ j_cap) ./ convert(T, sqrt(2)), std_axes, PE2D.Manifold(r_c1, PE2D.Axes(convert(T, -pi / 4)), (r_c1 - r_c1 / 2) .* (i_cap .+ j_cap) ./ convert(T, sqrt(2)))),

            # rotated_axes
            (c1, c2, origin, rotated_axes, PE2D.Manifold(r_c1 + r_c2, PE2D.rotate_minus_90(std_axes), (r_c1 - (r_c1 + r_c2) / 2) .* i_cap)),

            (c1, c2, r_c1 * i_cap, rotated_axes, PE2D.Manifold(r_c2, PE2D.rotate_minus_90(std_axes), (r_c1 - r_c2 / 2) .* i_cap)),
            (c1, c2, r_c2 * i_cap, rotated_axes, PE2D.Manifold(r_c1, PE2D.rotate_minus_90(std_axes), (r_c1 - r_c1 / 2) .* i_cap)),

            (c1, c2, r_c1 * j_cap, rotated_axes, PE2D.Manifold(r_c2, std_axes, (r_c1 - r_c2 / 2) .* j_cap)),
            (c1, c2, r_c2 * j_cap, rotated_axes, PE2D.Manifold(r_c1, std_axes, (r_c1 - r_c1 / 2) .* j_cap)),

            (c1, c2, r_c1 .* (i_cap .+ j_cap) ./ convert(T, sqrt(2)), rotated_axes, PE2D.Manifold(r_c2, PE2D.Axes(convert(T, -pi / 4)), (r_c1 - r_c2 / 2) .* (i_cap .+ j_cap) ./ convert(T, sqrt(2)))),
            (c1, c2, r_c2 .* (i_cap .+ j_cap) ./ convert(T, sqrt(2)), rotated_axes, PE2D.Manifold(r_c1, PE2D.Axes(convert(T, -pi / 4)), (r_c1 - r_c1 / 2) .* (i_cap .+ j_cap) ./ convert(T, sqrt(2)))),
            ]

            test_manifold_list(manifold_list)
        end

        @testset "Rect2D vs. Circle" begin
            manifold_list = [
            # std_axes
            (r1, c1, -(half_width_r1 + r_c1 - d) .* i_cap, std_axes, PE2D.Manifold(d, PE2D.rotate_90(std_axes), -(half_width_r1 + r_c1 - d) .* i_cap .+ (r_c1 - d / 2) .* i_cap)),
            (r1, c1, -(half_width_r1 + d) .* i_cap, std_axes, PE2D.Manifold(r_c1 - d, PE2D.rotate_90(std_axes), -(half_width_r1 + d) .* i_cap .+ (r_c1 - (r_c1 - d) / 2) .* i_cap)),
            (r1, c1, -(half_width_r1 - d) .* i_cap, std_axes, PE2D.Manifold(r_c1 + d, PE2D.rotate_90(std_axes), -(half_width_r1 - d) .* i_cap .+ (r_c1 - (r_c1 + d) / 2) .* i_cap)),
            (r1, c1, (half_width_r1 - d) .* i_cap, std_axes, PE2D.Manifold(r_c1 + d, PE2D.rotate_minus_90(std_axes), (half_width_r1 - d) .* i_cap .+ (r_c1 - (r_c1 + d) / 2) .* -i_cap)),
            (r1, c1, (half_width_r1 + d) .* i_cap, std_axes, PE2D.Manifold(r_c1 - d, PE2D.rotate_minus_90(std_axes), (half_width_r1 + d) .* i_cap .+ (r_c1 - (r_c1 - d) / 2) .* -i_cap)),
            (r1, c1, (half_width_r1 + r_c1 - d) .* i_cap, std_axes, PE2D.Manifold(d, PE2D.rotate_minus_90(std_axes), (half_width_r1 + r_c1 - d) .* i_cap .+ (r_c1 - d / 2) .* -i_cap)),

            (r1, c1, -(half_height_r1 + r_c1 - d) .* j_cap, std_axes, PE2D.Manifold(d, PE2D.rotate_180(std_axes), -(half_height_r1 + r_c1 - d) .* j_cap .+ (r_c1 - d / 2) .* j_cap)),
            (r1, c1, -(half_height_r1 + d) .* j_cap, std_axes, PE2D.Manifold(r_c1 - d, PE2D.rotate_180(std_axes), -(half_height_r1 + d) .* j_cap .+ (r_c1 - (r_c1 - d) / 2) .* j_cap)),
            (r1, c1, -(half_height_r1 - d) .* j_cap, std_axes, PE2D.Manifold(r_c1 + d, PE2D.rotate_180(std_axes), -(half_height_r1 - d) .* j_cap .+ (r_c1 - (r_c1 + d) / 2) .* j_cap)),
            (r1, c1, -d .* j_cap, std_axes, PE2D.Manifold(half_height_r1 + r_c1 - d, PE2D.rotate_180(std_axes), -d .* j_cap .+ (r_c1 - (half_height_r1 + r_c1 - d) / 2) .* j_cap)),
            (r1, c1, d .* j_cap, std_axes, PE2D.Manifold(half_height_r1 + r_c1 - d, std_axes, d .* j_cap .+ (r_c1 - (half_height_r1 + r_c1 - d) / 2) .* -j_cap)),
            (r1, c1, (half_height_r1 - d) .* j_cap, std_axes, PE2D.Manifold(r_c1 + d, std_axes, (half_height_r1 - d) .* j_cap .+ (r_c1 - (r_c1 + d) / 2) .* -j_cap)),
            (r1, c1, (half_height_r1 + d) .* j_cap, std_axes, PE2D.Manifold(r_c1 - d, std_axes, (half_height_r1 + d) .* j_cap .+ (r_c1 - (r_c1 - d) / 2) .* -j_cap)),
            (r1, c1, (half_height_r1 + r_c1 - d) .* j_cap, std_axes, PE2D.Manifold(d, std_axes, (half_height_r1 + r_c1 - d) .* j_cap .+ (r_c1 - d / 2) .* -j_cap)),

            (r1, c1, half_width_r1/2 .* i_cap .+ (half_height_r1 + r_c1 - d) .* j_cap, std_axes, PE2D.Manifold(d, std_axes, half_width_r1/2 .* i_cap .+ (half_height_r1 + r_c1 - d) .* j_cap .+ (r_c1 - d / 2) .* -j_cap)),

            # rotated_axes
            (r1, c1, -(half_width_r1 + r_c1 - d) .* i_cap, rotated_axes, PE2D.Manifold(d, PE2D.rotate_90(std_axes), -(half_width_r1 + r_c1 - d) .* i_cap .+ (r_c1 - d / 2) .* i_cap)),
            (r1, c1, -(half_width_r1 + d) .* i_cap, rotated_axes, PE2D.Manifold(r_c1 - d, PE2D.rotate_90(std_axes), -(half_width_r1 + d) .* i_cap .+ (r_c1 - (r_c1 - d) / 2) .* i_cap)),
            (r1, c1, -(half_width_r1 - d) .* i_cap, rotated_axes, PE2D.Manifold(r_c1 + d, PE2D.rotate_90(std_axes), -(half_width_r1 - d) .* i_cap .+ (r_c1 - (r_c1 + d) / 2) .* i_cap)),
            (r1, c1, (half_width_r1 - d) .* i_cap, rotated_axes, PE2D.Manifold(r_c1 + d, PE2D.rotate_minus_90(std_axes), (half_width_r1 - d) .* i_cap .+ (r_c1 - (r_c1 + d) / 2) .* -i_cap)),
            (r1, c1, (half_width_r1 + d) .* i_cap, rotated_axes, PE2D.Manifold(r_c1 - d, PE2D.rotate_minus_90(std_axes), (half_width_r1 + d) .* i_cap .+ (r_c1 - (r_c1 - d) / 2) .* -i_cap)),
            (r1, c1, (half_width_r1 + r_c1 - d) .* i_cap, rotated_axes, PE2D.Manifold(d, PE2D.rotate_minus_90(std_axes), (half_width_r1 + r_c1 - d) .* i_cap .+ (r_c1 - d / 2) .* -i_cap)),

            (r1, c1, -(half_height_r1 + r_c1 - d) .* j_cap, rotated_axes, PE2D.Manifold(d, PE2D.rotate_180(std_axes), -(half_height_r1 + r_c1 - d) .* j_cap .+ (r_c1 - d / 2) .* j_cap)),
            (r1, c1, -(half_height_r1 + d) .* j_cap, rotated_axes, PE2D.Manifold(r_c1 - d, PE2D.rotate_180(std_axes), -(half_height_r1 + d) .* j_cap .+ (r_c1 - (r_c1 - d) / 2) .* j_cap)),
            (r1, c1, -(half_height_r1 - d) .* j_cap, rotated_axes, PE2D.Manifold(r_c1 + d, PE2D.rotate_180(std_axes), -(half_height_r1 - d) .* j_cap .+ (r_c1 - (r_c1 + d) / 2) .* j_cap)),
            (r1, c1, -d .* j_cap, rotated_axes, PE2D.Manifold(half_height_r1 + r_c1 - d, PE2D.rotate_180(std_axes), -d .* j_cap .+ (r_c1 - (half_height_r1 + r_c1 - d) / 2) .* j_cap)),
            (r1, c1, d .* j_cap, rotated_axes, PE2D.Manifold(half_height_r1 + r_c1 - d, std_axes, d .* j_cap .+ (r_c1 - (half_height_r1 + r_c1 - d) / 2) .* -j_cap)),
            (r1, c1, (half_height_r1 - d) .* j_cap, rotated_axes, PE2D.Manifold(r_c1 + d, std_axes, (half_height_r1 - d) .* j_cap .+ (r_c1 - (r_c1 + d) / 2) .* -j_cap)),
            (r1, c1, (half_height_r1 + d) .* j_cap, rotated_axes, PE2D.Manifold(r_c1 - d, std_axes, (half_height_r1 + d) .* j_cap .+ (r_c1 - (r_c1 - d) / 2) .* -j_cap)),
            (r1, c1, (half_height_r1 + r_c1 - d) .* j_cap, rotated_axes, PE2D.Manifold(d, std_axes, (half_height_r1 + r_c1 - d) .* j_cap .+ (r_c1 - d / 2) .* -j_cap)),

            (r1, c1, half_width_r1/2 .* i_cap .+ (half_height_r1 + r_c1 - d) .* j_cap, rotated_axes, PE2D.Manifold(d, std_axes, half_width_r1/2 .* i_cap .+ (half_height_r1 + r_c1 - d) .* j_cap .+ (r_c1 - d / 2) .* -j_cap)),
            ]

            test_manifold_list(manifold_list)
        end

        # @testset "Rect2D vs. Rect2D" begin
            # (r1, r2, origin, std_axes, PE2D.Manifold(half_height_r1 + half_height_r2, PE2D.rotate_minus_90(std_axes), (r_c1 - (r_c1 + r_c2) / 2) .* i_cap)),
            # manifold_list = [(GB.HyperRectangle(1.0f0, 2.0f0, 3.0f0, 4.0f0), GB.HyperRectangle(0.0f0, 0.0f0, 1.1f0, 5.0f0), PE2D.Manifold(0.1f0, GB.Vec(-1.0f0, 0.0f0), GB.Vec(1.05f0, 3.5f0))),
                             # (GB.HyperRectangle(1.0f0, 2.0f0, 3.0f0, 4.0f0), GB.HyperRectangle(2.0f0, 0.0f0, 5.0f0, 2.1f0), PE2D.Manifold(0.1f0, GB.Vec(0.0f0, -1.0f0), GB.Vec(3.0f0, 2.05f0))),
                             # (GB.HyperRectangle(1.0f0, 2.0f0, 3.0f0, 4.0f0), GB.HyperRectangle(2.0f0, 0.0f0, 5.0f0, 2.0f0), PE2D.Manifold(0.0f0, GB.Vec(0.0f0, -1.0f0), GB.Vec(3.0f0, 2.0f0)))]
            # test_manifold_list(manifold_list)
        # end
    end

    # @testset "Simulation: frictionless collision of rotating circles" begin
        # T = Float32
        # NUM_ITER = 500
        # DT = 1 / 60
        # FRAME_RATE = 1 / DT

        # shape1 = GB.HyperSphere(GB.Point2{T}(5, 1), one(T))
        # material_data1 = PE2D.MaterialData{T}()
        # mass_data1 = PE2D.MassData(material_data1.density, shape1)
        # position_accumulator1 = PE2D.Accumulator(GB.Vec2{T}(5, 1), zero(GB.Vec2{T}))
        # velocity_accumulator1 = PE2D.Accumulator(GB.Vec2{T}(0, 1), zero(GB.Vec2{T}))
        # force_accumulator1 = PE2D.Accumulator(zero(GB.Vec2{T}), zero(GB.Vec2{T}))
        # inertia_data1 = PE2D.InertiaData(material_data1.density, shape1)
        # angle1 = zero(T)
        # angle_accumulator1 = PE2D.Accumulator(angle1, zero(T))
        # angular_velocity_accumulator1 = PE2D.Accumulator(one(T), zero(T))
        # torque_accumulator1 = PE2D.Accumulator(zero(T), zero(T))
        # direction1 = GB.Vec2{T}(cos(angle1), sin(angle1))
        # body1 = PE2D.RigidBody(shape1, material_data1, mass_data1, position_accumulator1, velocity_accumulator1, force_accumulator1, inertia_data1, angle_accumulator1, angular_velocity_accumulator1, torque_accumulator1, direction1)

        # shape2 = GB.HyperSphere(GB.Point2{T}(1, 5), one(T))
        # material_data2 = PE2D.MaterialData{T}()
        # mass_data2 = PE2D.MassData(material_data2.density, shape2)
        # position_accumulator2 = PE2D.Accumulator(GB.Vec2{T}(1, 5), zero(GB.Vec2{T}))
        # velocity_accumulator2 = PE2D.Accumulator(GB.Vec2{T}(1, 0), zero(GB.Vec2{T}))
        # force_accumulator2 = PE2D.Accumulator(zero(GB.Vec2{T}), zero(GB.Vec2{T}))
        # inertia_data2 = PE2D.InertiaData(material_data2.density, shape2)
        # angle2 = zero(T)
        # angle_accumulator2 = PE2D.Accumulator(angle2, zero(T))
        # angular_velocity_accumulator2 = PE2D.Accumulator(one(T), zero(T))
        # torque_accumulator2 = PE2D.Accumulator(zero(T), zero(T))
        # direction2 = GB.Vec2{T}(cos(angle2), sin(angle2))
        # body2 = PE2D.RigidBody(shape2, material_data2, mass_data2, position_accumulator2, velocity_accumulator2, force_accumulator2, inertia_data2, angle_accumulator2, angular_velocity_accumulator2, torque_accumulator2, direction2)

        # bodies = [body1, body2]
        # world = PE2D.World(bodies)
        # PE2D.simulate!(world, NUM_ITER, DT)
    # end

end
