import PhysicsEngine2D
import PhysicsEngine2D: PE2D
import GeometryBasics
const GB = GeometryBasics
import LinearAlgebra
const LA = LinearAlgebra
using Test

function test_collision_list(collision_list)
    for (a, b, pos_ba, axes_ba, value) in collision_list
        # @show a
        # @show b
        # @show pos_ba
        # @show axes_ba
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
    theta_45 = convert(T, π / 4)
    unit_45 = (i_cap .+ j_cap) / convert(T, sqrt(2))

    l1 = GB.Line(convert(GB.Point, -i_cap), convert(GB.Point, i_cap))
    p1_l1 = PE2D.get_point(l1, 1)
    p2_l1 = PE2D.get_point(l1, 2)
    half_width_l1 = p2_l1[1]

    l2 = GB.Line(convert(GB.Point, -2 .* i_cap), convert(GB.Point, 2 .* i_cap))
    p1_l2 = PE2D.get_point(l2, 1)
    p2_l2 = PE2D.get_point(l2, 2)
    half_width_l2 = p2_l2[1]

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
        @testset "Line segment vs. Line segment" begin
            collision_list = [
            # std_axes
            (l1, l2, origin, std_axes, true),

            (l1, l2, (half_width_l1 + half_width_l2 + d) .* -i_cap, std_axes, false),
            (l1, l2, (half_width_l1 + half_width_l2 - d) .* -i_cap, std_axes, true),
            (l1, l2, (half_width_l1 + half_width_l2 - d) .* i_cap, std_axes, true),
            (l1, l2, (half_width_l1 + half_width_l2 + d) .* i_cap, std_axes, false),

            (l1, l2, d .* -j_cap, std_axes, false),
            (l1, l2, d .* j_cap, std_axes, false),

            # rotated_axes
            (l1, l2, origin, rotated_axes, true),

            (l1, l2, (half_width_l1 + half_width_l2 + d) .* -i_cap, rotated_axes, false),
            (l1, l2, (half_width_l1 + half_width_l2 - d) .* -i_cap, rotated_axes, false),
            (l1, l2, (half_width_l1 + d) .* -i_cap, rotated_axes, false),
            (l1, l2, (half_width_l1 - d) .* -i_cap, rotated_axes, true),
            (l1, l2, (half_width_l1 - d) .* i_cap, rotated_axes, true),
            (l1, l2, (half_width_l1 + d) .* i_cap, rotated_axes, false),
            (l1, l2, (half_width_l1 + half_width_l2 - d) .* i_cap, rotated_axes, false),
            (l1, l2, (half_width_l1 + half_width_l2 + d) .* i_cap, rotated_axes, false),

            (l2, l1, (half_width_l1 * sin(theta) + d) .* -j_cap, rotated_axes, false),
            (l2, l1, (half_width_l1 * sin(theta) - d) .* -j_cap, rotated_axes, true),
            (l2, l1, d .* -j_cap, rotated_axes, true),
            (l2, l1, d .* j_cap, rotated_axes, true),
            (l2, l1, (half_width_l1 * sin(theta) - d) .* j_cap, rotated_axes, true),
            (l2, l1, (half_width_l1 * sin(theta) + d) .* j_cap, rotated_axes, false),
            ]

            test_collision_list(collision_list)
        end

        @testset "Circle vs. Point2" begin
            collision_list = [
            # std_axes
            (c1, origin, origin, std_axes, true),

            (c1, origin, (r_c1 + d) .* -i_cap, std_axes, false),
            (c1, origin, (r_c1 - d) .* -i_cap, std_axes, true),
            (c1, origin, (r_c1 - d) .* i_cap, std_axes, true),
            (c1, origin, (r_c1 + d) .* i_cap, std_axes, false),

            (c1, origin, (r_c1 + d) .* -j_cap, std_axes, false),
            (c1, origin, (r_c1 - d) .* -j_cap, std_axes, true),
            (c1, origin, (r_c1 - d) .* j_cap, std_axes, true),
            (c1, origin, (r_c1 + d) .* j_cap, std_axes, false),

            (c1, origin, (r_c1 + d) .* unit_45, std_axes, false),
            (c1, origin, (r_c1 - d) .* unit_45, std_axes, true),

            # reverse check with std_axes
            (origin, c1, origin, std_axes, true),

            (origin, c1, (r_c1 + d) .* -i_cap, std_axes, false),
            (origin, c1, (r_c1 - d) .* -i_cap, std_axes, true),
            (origin, c1, (r_c1 - d) .* i_cap, std_axes, true),
            (origin, c1, (r_c1 + d) .* i_cap, std_axes, false),

            (origin, c1, (r_c1 + d) .* -j_cap, std_axes, false),
            (origin, c1, (r_c1 - d) .* -j_cap, std_axes, true),
            (origin, c1, (r_c1 - d) .* j_cap, std_axes, true),
            (origin, c1, (r_c1 + d) .* j_cap, std_axes, false),

            (origin, c1, (r_c1 + d) .* unit_45, std_axes, false),
            (origin, c1, (r_c1 - d) .* unit_45, std_axes, true),
            ]

            test_collision_list(collision_list)
        end

        @testset "Circle vs. Line segment" begin
            collision_list = [
            # std_axes
            (l1, c2, origin, std_axes, true),
            (l2, c1, origin, std_axes, true),

            (l2, c1, (half_width_l2 + r_c1 + d) .* -i_cap, std_axes, false),
            (l2, c1, (half_width_l2 + r_c1 - d) .* -i_cap, std_axes, true),
            (l2, c1, (half_width_l2 + r_c1 - d) .* i_cap, std_axes, true),
            (l2, c1, (half_width_l2 + r_c1 + d) .* i_cap, std_axes, false),

            (l2, c1, (r_c1 + d) .* -j_cap, std_axes, false),
            (l2, c1, (r_c1 - d) .* -j_cap, std_axes, true),
            (l2, c1, (r_c1 - d) .* j_cap, std_axes, true),
            (l2, c1, (r_c1 + d) .* j_cap, std_axes, false),

            # reverse check with std_axes
            (c2, l1, origin, std_axes, true),
            (c1, l2, origin, std_axes, true),

            (c1, l2, (half_width_l2 + r_c1 + d) .* -i_cap, std_axes, false),
            (c1, l2, (half_width_l2 + r_c1 - d) .* -i_cap, std_axes, true),
            (c1, l2, (half_width_l2 + r_c1 - d) .* i_cap, std_axes, true),
            (c1, l2, (half_width_l2 + r_c1 + d) .* i_cap, std_axes, false),

            (c1, l2, (r_c1 + d) .* -j_cap, std_axes, false),
            (c1, l2, (r_c1 - d) .* -j_cap, std_axes, true),
            (c1, l2, (r_c1 - d) .* j_cap, std_axes, true),
            (c1, l2, (r_c1 + d) .* j_cap, std_axes, false),

            # reverse check with rotated_axes
            (c2, l1, origin, rotated_axes, true),
            (c1, l2, origin, rotated_axes, true),

            (c1, l1, (sqrt(r_c1 ^ 2 - (half_width_l1 * sin(theta)) ^ 2) + half_width_l1 * cos(theta) + d) .* -i_cap, rotated_axes, false),
            (c1, l1, (sqrt(r_c1 ^ 2 - (half_width_l1 * sin(theta)) ^ 2) + half_width_l1 * cos(theta) - d) .* -i_cap, rotated_axes, true),
            (c1, l1, (sqrt(r_c1 ^ 2 - (half_width_l1 * sin(theta)) ^ 2) + half_width_l1 * cos(theta) - d) .* i_cap, rotated_axes, true),
            (c1, l1, (sqrt(r_c1 ^ 2 - (half_width_l1 * sin(theta)) ^ 2) + half_width_l1 * cos(theta) + d) .* i_cap, rotated_axes, false),

            (c1, l2, (r_c1 / cos(theta) + d) .* -j_cap, rotated_axes, false),
            (c1, l2, (r_c1 / cos(theta) - d) .* -j_cap, rotated_axes, true),
            (c1, l2, (r_c1 / cos(theta) - d) .* j_cap, rotated_axes, true),
            (c1, l2, (r_c1 / cos(theta) + d) .* j_cap, rotated_axes, false),
            ]

            test_collision_list(collision_list)
        end

        @testset "Circle vs. Circle" begin
            collision_list = [
            # std_axes
            (c1, c2, origin, std_axes, true),

            (c1, c2, (r_c1 + r_c2 + d) .* -i_cap, std_axes, false),
            (c1, c2, (r_c1 + r_c2 - d) .* -i_cap, std_axes, true),
            (c1, c2, (r_c1 + r_c2 - d) .* i_cap, std_axes, true),
            (c1, c2, (r_c1 + r_c2 + d) .* i_cap, std_axes, false),

            (c1, c2, (r_c1 + r_c2 + d) .* -j_cap, std_axes, false),
            (c1, c2, (r_c1 + r_c2 - d) .* -j_cap, std_axes, true),
            (c1, c2, (r_c1 + r_c2 - d) .* j_cap, std_axes, true),
            (c1, c2, (r_c1 + r_c2 + d) .* j_cap, std_axes, false),

            (c1, c2, (r_c1 + r_c2 + d) .* unit_45, std_axes, false),
            (c1, c2, (r_c1 + r_c2 - d) .* unit_45, std_axes, true),
            ]

            test_collision_list(collision_list)
        end

        @testset "Rect2D vs. Point2" begin
            collision_list = [
            # std_axes
            (r1, origin, origin, std_axes, true),

            (r1, origin, (half_width_r1 + d) .* -i_cap, std_axes, false),
            (r1, origin, (half_width_r1 - d) .* -i_cap, std_axes, true),
            (r1, origin, (half_width_r1 - d) .* i_cap, std_axes, true),
            (r1, origin, (half_width_r1 + d) .* i_cap, std_axes, false),

            (r1, origin, (half_height_r1 + d) .* -j_cap, std_axes, false),
            (r1, origin, (half_height_r1 - d) .* -j_cap, std_axes, true),
            (r1, origin, (half_height_r1 - d) .* j_cap, std_axes, true),
            (r1, origin, (half_height_r1 + d) .* j_cap, std_axes, false),

            (r1, origin, top_right_r1 .+ d, std_axes, false),
            (r1, origin, top_right_r1 .- d, std_axes, true),

            # reverse check with std_axes
            (origin, r1, origin, std_axes, true),

            (origin, r1, (half_width_r1 + d) .* -i_cap, std_axes, false),
            (origin, r1, (half_width_r1 - d) .* -i_cap, std_axes, true),
            (origin, r1, (half_width_r1 - d) .* i_cap, std_axes, true),
            (origin, r1, (half_width_r1 + d) .* i_cap, std_axes, false),

            (origin, r1, (half_height_r1 + d) .* -j_cap, std_axes, false),
            (origin, r1, (half_height_r1 - d) .* -j_cap, std_axes, true),
            (origin, r1, (half_height_r1 - d) .* j_cap, std_axes, true),
            (origin, r1, (half_height_r1 + d) .* j_cap, std_axes, false),

            (origin, r1, top_right_r1 .+ d, std_axes, false),
            (origin, r1, top_right_r1 .- d, std_axes, true),

            # reverse check with rotated_axes
            (origin, r1, origin, rotated_axes, true),

            (origin, r1, (half_height_r1 / sin(theta) + d) .* -i_cap, rotated_axes, false),
            (origin, r1, (half_height_r1 / sin(theta) - d) .* -i_cap, rotated_axes, true),
            (origin, r1, (half_height_r1 / sin(theta) - d) .* i_cap, rotated_axes, true),
            (origin, r1, (half_height_r1 / sin(theta) + d) .* i_cap, rotated_axes, false),

            (origin, r1, (half_height_r1 / cos(theta) + d) .* -j_cap, rotated_axes, false),
            (origin, r1, (half_height_r1 / cos(theta) - d) .* -j_cap, rotated_axes, true),
            (origin, r1, (half_height_r1 / cos(theta) - d) .* j_cap, rotated_axes, true),
            (origin, r1, (half_height_r1 / cos(theta) + d) .* j_cap, rotated_axes, false),
            ]

            test_collision_list(collision_list)
        end

        @testset "Rect2D vs. Line segment" begin
            collision_list = [
            # std_axes
            (r1, l1, origin, std_axes, true),
            (r1, l2, origin, std_axes, true),

            (r1, l1, (half_width_r1 + half_width_l1 + d) .* -i_cap, std_axes, false),
            (r1, l1, (half_width_r1 + half_width_l1 - d) .* -i_cap, std_axes, true),
            (r1, l1, (half_width_r1 + half_width_l1 - d) .* i_cap, std_axes, true),
            (r1, l1, (half_width_r1 + half_width_l1 + d) .* i_cap, std_axes, false),

            (r1, l1, (half_height_r1 + d) .* -j_cap, std_axes, false),
            (r1, l1, (half_height_r1 - d) .* -j_cap, std_axes, true),
            (r1, l1, (half_height_r1 - d) .* j_cap, std_axes, true),
            (r1, l1, (half_height_r1 + d) .* j_cap, std_axes, false),

            # rotated_axes
            (r1, l1, origin, rotated_axes, true),
            (r1, l2, origin, rotated_axes, true),

            (r2, l1, (half_width_r2 + half_width_l1 * cos(theta) + d) .* -i_cap, rotated_axes, false),
            (r2, l1, (half_width_r2 + half_width_l1 * cos(theta) - d) .* -i_cap, rotated_axes, true),
            (r2, l1, (half_width_r2 + half_width_l1 * cos(theta) - d) .* i_cap, rotated_axes, true),
            (r2, l1, (half_width_r2 + half_width_l1 * cos(theta) + d) .* i_cap, rotated_axes, false),

            (r2, l1, (half_height_r2 + half_width_l1 * sin(theta) + d) .* -j_cap, rotated_axes, false),
            (r2, l1, (half_height_r2 + half_width_l1 * sin(theta) - d) .* -j_cap, rotated_axes, true),
            (r2, l1, (half_height_r2 + half_width_l1 * sin(theta) - d) .* j_cap, rotated_axes, true),
            (r2, l1, (half_height_r2 + half_width_l1 * sin(theta) + d) .* j_cap, rotated_axes, false),

            # reverse check with std_axes
            (l1, r1, origin, std_axes, true),
            (l2, r1, origin, std_axes, true),

            (l1, r1, (half_width_r1 + half_width_l1 + d) .* -i_cap, std_axes, false),
            (l1, r1, (half_width_r1 + half_width_l1 - d) .* -i_cap, std_axes, true),
            (l1, r1, (half_width_r1 + half_width_l1 - d) .* i_cap, std_axes, true),
            (l1, r1, (half_width_r1 + half_width_l1 + d) .* i_cap, std_axes, false),

            (l1, r1, (half_height_r1 + d) .* -j_cap, std_axes, false),
            (l1, r1, (half_height_r1 - d) .* -j_cap, std_axes, true),
            (l1, r1, (half_height_r1 - d) .* j_cap, std_axes, true),
            (l1, r1, (half_height_r1 + d) .* j_cap, std_axes, false),

            # reverse check with rotated_axes
            (l1, r1, origin, rotated_axes, true),
            (l2, r1, origin, rotated_axes, true),

            (l2, r2, (half_height_r2 / sin(theta) + half_width_l2 + d) .* -i_cap, rotated_axes, false),
            (l2, r2, (half_height_r2 / sin(theta) + half_width_l2 - d) .* -i_cap, rotated_axes, true),
            (l2, r2, (half_height_r2 / sin(theta) + half_width_l2 - d) .* i_cap, rotated_axes, true),
            (l2, r2, (half_height_r2 / sin(theta) + half_width_l2 + d) .* i_cap, rotated_axes, false),

            (l2, r2, (half_width_r2 * sin(theta) + half_height_r2 * cos(theta) + d) .* -j_cap, rotated_axes, false),
            (l2, r2, (half_width_r2 * sin(theta) + half_height_r2 * cos(theta) - d) .* -j_cap, rotated_axes, true),
            (l2, r2, (half_width_r2 * sin(theta) + half_height_r2 * cos(theta) - d) .* j_cap, rotated_axes, true),
            (l2, r2, (half_width_r2 * sin(theta) + half_height_r2 * cos(theta) + d) .* j_cap, rotated_axes, false),
            ]

            test_collision_list(collision_list)
        end

        @testset "Rect2D vs. Circle" begin
            top_right = PE2D.get_top_right(r1)
            theta_r1 = atan(top_right[2], top_right[1])

            collision_list = [
            # std_axes
            (r1, c1, origin, std_axes, true),

            (r1, c1, (half_width_r1 + r_c1 + d) .* -i_cap, std_axes, false),
            (r1, c1, (half_width_r1 + r_c1 - d) .* -i_cap, std_axes, true),
            (r1, c1, (half_width_r1 + r_c1 - d) .* i_cap, std_axes, true),
            (r1, c1, (half_width_r1 + r_c1 + d) .* i_cap, std_axes, false),

            (r1, c1, (half_height_r1 + r_c1 + d) .* -j_cap, std_axes, false),
            (r1, c1, (half_height_r1 + r_c1 - d) .* -j_cap, std_axes, true),
            (r1, c1, (half_height_r1 + r_c1 - d) .* j_cap, std_axes, true),
            (r1, c1, (half_height_r1 + r_c1 + d) .* -j_cap, std_axes, false),

            (r1, c1, top_right_r1 .+ (r_c1 + d) .* unit_45, std_axes, false),
            (r1, c1, top_right_r1 .+ (r_c1 - d) .* unit_45, std_axes, true),

            # reverse check with std_axes
            (c1, r1, origin, std_axes, true),

            (c1, r1, (half_width_r1 + r_c1 + d) .* -i_cap, std_axes, false),
            (c1, r1, (half_width_r1 + r_c1 - d) .* -i_cap, std_axes, true),
            (c1, r1, (half_width_r1 + r_c1 - d) .* i_cap, std_axes, true),
            (c1, r1, (half_width_r1 + r_c1 + d) .* i_cap, std_axes, false),

            (c1, r1, (half_height_r1 + r_c1 + d) .* -j_cap, std_axes, false),
            (c1, r1, (half_height_r1 + r_c1 - d) .* -j_cap, std_axes, true),
            (c1, r1, (half_height_r1 + r_c1 - d) .* j_cap, std_axes, true),
            (c1, r1, (half_height_r1 + r_c1 + d) .* -j_cap, std_axes, false),

            (c1, r1, top_right_r1 .+ (r_c1 + d) .* unit_45, std_axes, false),
            (c1, r1, top_right_r1 .+ (r_c1 - d) .* unit_45, std_axes, true),

            # reverse check with rotated_axes
            (c2, r2, origin, rotated_axes, true),

            (c2, r2, (sqrt(r_c2 ^ 2 - (half_width_r2 * sin(theta) - half_height_r2 * cos(theta)) ^ 2) + half_width_r2 * cos(theta) + half_height_r2 * sin(theta) + d) .* -i_cap, rotated_axes, false),
            (c2, r2, (sqrt(r_c2 ^ 2 - (half_width_r2 * sin(theta) - half_height_r2 * cos(theta)) ^ 2) + half_width_r2 * cos(theta) + half_height_r2 * sin(theta) - d) .* -i_cap, rotated_axes, true),
            (c2, r2, (sqrt(r_c2 ^ 2 - (half_width_r2 * sin(theta) - half_height_r2 * cos(theta)) ^ 2) + half_width_r2 * cos(theta) + half_height_r2 * sin(theta) - d) .* i_cap, rotated_axes, true),
            (c2, r2, (sqrt(r_c2 ^ 2 - (half_width_r2 * sin(theta) - half_height_r2 * cos(theta)) ^ 2) + half_width_r2 * cos(theta) + half_height_r2 * sin(theta) + d) .* i_cap, rotated_axes, false),

            (c2, r2, ((r_c2 + half_height_r2) / cos(theta) + d) .* -j_cap, rotated_axes, false),
            (c2, r2, ((r_c2 + half_height_r2) / cos(theta) - d) .* -j_cap, rotated_axes, true),
            (c2, r2, ((r_c2 + half_height_r2) / cos(theta) - d) .* j_cap, rotated_axes, true),
            (c2, r2, ((r_c2 + half_height_r2) / cos(theta) + d) .* j_cap, rotated_axes, false),
            ]

            test_collision_list(collision_list)
        end

        @testset "Rect2D vs. Rect2D" begin
            collision_list = [
            # std_axes
            (r2, r1, origin, std_axes, true),

            (r2, r1, (half_width_r1 + half_width_r2 + d) .* -i_cap, std_axes, false),
            (r2, r1, (half_width_r1 + half_width_r2 - d) .* -i_cap, std_axes, true),
            (r2, r1, (half_width_r1 + half_width_r2 - d) .* i_cap, std_axes, true),
            (r2, r1, (half_width_r1 + half_width_r2 + d) .* i_cap, std_axes, false),

            (r2, r1, (half_height_r1 + half_height_r2 + d) .* -j_cap, std_axes, false),
            (r2, r1, (half_height_r1 + half_height_r2 - d) .* -j_cap, std_axes, true),
            (r2, r1, (half_height_r1 + half_height_r2 - d) .* j_cap, std_axes, true),
            (r2, r1, (half_height_r1 + half_height_r2 + d) .* j_cap, std_axes, false),

            (r2, r1, top_right_r1 .+ top_right_r2 .+ d, std_axes, false),
            (r2, r1, top_right_r1 .+ top_right_r2 .- d, std_axes, true),

            # rotated_axes
            (r2, r1, origin, rotated_axes, true),

            (r2, r1, (half_width_r2 + half_width_r1 * cos(theta) + half_height_r1 * sin(theta) + d) .* -i_cap, rotated_axes, false),
            (r2, r1, (half_width_r2 + half_width_r1 * cos(theta) + half_height_r1 * sin(theta) - d) .* -i_cap, rotated_axes, true),
            (r2, r1, (half_width_r2 + half_width_r1 * cos(theta) + half_height_r1 * sin(theta) - d) .* i_cap, rotated_axes, true),
            (r2, r1, (half_width_r2 + half_width_r1 * cos(theta) + half_height_r1 * sin(theta) + d) .* i_cap, rotated_axes, false),

            (r2, r1, (half_height_r2 + half_width_r1 * sin(theta) + half_height_r1 * cos(theta) + d) .* -j_cap, rotated_axes, false),
            (r2, r1, (half_height_r2 + half_width_r1 * sin(theta) + half_height_r1 * cos(theta) - d) .* -j_cap, rotated_axes, true),
            (r2, r1, (half_height_r2 + half_width_r1 * sin(theta) + half_height_r1 * cos(theta) - d) .* j_cap, rotated_axes, true),
            (r2, r1, (half_height_r2 + half_width_r1 * sin(theta) + half_height_r1 * cos(theta) + d) .* j_cap, rotated_axes, false),
            ]

            test_collision_list(collision_list)
        end
    end

    @testset "Manifold generation" begin
        @testset "Circle vs. Circle" begin
            manifold_list = [
            # std_axes
            (c1, c2, (r_c1 + r_c2 - d) .* -i_cap, std_axes, PE2D.Manifold(d, PE2D.rotate_90(std_axes), (r_c1 - d / 2) .* -i_cap)),
            (c1, c2, r_c2 .* -i_cap, std_axes, PE2D.Manifold(r_c1, PE2D.rotate_90(std_axes), (r_c1 / 2) .* -i_cap)),
            (c1, c2, r_c2 .* i_cap, std_axes, PE2D.Manifold(r_c1, PE2D.rotate_minus_90(std_axes), (r_c1 / 2) .* i_cap)),
            (c1, c2, (r_c1 + r_c2 - d) .* i_cap, std_axes, PE2D.Manifold(d, PE2D.rotate_minus_90(std_axes), (r_c1 - d / 2) .* i_cap)),

            (c1, c2, (r_c1 + r_c2 - d) .* -j_cap, std_axes, PE2D.Manifold(d, PE2D.rotate_180(std_axes), (r_c1 - d / 2) .* -j_cap)),
            (c1, c2, r_c2 .* -j_cap, std_axes, PE2D.Manifold(r_c1, PE2D.rotate_180(std_axes), (r_c1 / 2) .* -j_cap)),
            (c1, c2, r_c2 .* j_cap, std_axes, PE2D.Manifold(r_c1, std_axes, (r_c1 / 2) .* j_cap)),
            (c1, c2, (r_c1 + r_c2 - d) .* j_cap, std_axes, PE2D.Manifold(d, std_axes, (r_c1 - d / 2) .* j_cap)),

            (c1, c2, (r_c1 + r_c2 - d) .* unit_45, std_axes, PE2D.Manifold(d, PE2D.Axes(-theta_45), (r_c1 - d / 2) .* unit_45)),
            (c1, c2, r_c2 .* unit_45, std_axes, PE2D.Manifold(r_c1, PE2D.Axes(-theta_45), (r_c1 / 2) .* unit_45)),
            ]

            test_manifold_list(manifold_list)
        end

        @testset "Rect2D vs. Circle" begin
            manifold_list = [
            # std_axes
            (r1, c1, (half_width_r1 + r_c1 - d) .* -i_cap, std_axes, PE2D.Manifold(d, PE2D.rotate_90(std_axes), (half_width_r1 - d / 2) .* -i_cap)),
            (r1, c1, (half_width_r1 + d) .* -i_cap, std_axes, PE2D.Manifold(r_c1 - d, PE2D.rotate_90(std_axes), (half_width_r1 - (r_c1 - d) / 2) .* -i_cap)),
            (r1, c1, (half_width_r1 - d) .* -i_cap, std_axes, PE2D.Manifold(r_c1 + d, PE2D.rotate_90(std_axes), (half_width_r1 - (r_c1 + d) / 2) .* -i_cap)),
            (r1, c1, (half_width_r1 - d) .* i_cap, std_axes, PE2D.Manifold(r_c1 + d, PE2D.rotate_minus_90(std_axes), (half_width_r1 - (r_c1 + d) / 2) .* i_cap)),
            (r1, c1, (half_width_r1 + d) .* i_cap, std_axes, PE2D.Manifold(r_c1 - d, PE2D.rotate_minus_90(std_axes), (half_width_r1 - (r_c1 - d) / 2) .* i_cap)),
            (r1, c1, (half_width_r1 + r_c1 - d) .* i_cap, std_axes, PE2D.Manifold(d, PE2D.rotate_minus_90(std_axes), (half_width_r1 - d / 2) .* i_cap)),

            (r1, c1, (half_height_r1 + r_c1 - d) .* -j_cap, std_axes, PE2D.Manifold(d, PE2D.rotate_180(std_axes), (half_height_r1 - d / 2) .* -j_cap)),
            (r1, c1, (half_height_r1 + d) .* -j_cap, std_axes, PE2D.Manifold(r_c1 - d, PE2D.rotate_180(std_axes), (half_height_r1 - (r_c1 - d) / 2) .* -j_cap)),
            (r1, c1, (half_height_r1 - d) .* -j_cap, std_axes, PE2D.Manifold(r_c1 + d, PE2D.rotate_180(std_axes), (half_height_r1 - (r_c1 + d) / 2) .* -j_cap)),
            (r1, c1, (half_height_r1 - d) .* j_cap, std_axes, PE2D.Manifold(r_c1 + d, std_axes, (half_height_r1 - (r_c1 + d) / 2) .* j_cap)),
            (r1, c1, (half_height_r1 + d) .* j_cap, std_axes, PE2D.Manifold(r_c1 - d, std_axes, (half_height_r1 - (r_c1 - d) / 2) .* j_cap)),
            (r1, c1, (half_height_r1 + r_c1 - d) .* j_cap, std_axes, PE2D.Manifold(d, std_axes, (half_height_r1 - d / 2) .* j_cap)),

            (r1, c1, top_right_r1 .+ (r_c1 - d) .* unit_45, std_axes, PE2D.Manifold(d, PE2D.Axes(-theta_45), top_right_r1 .+ (d / 2) .* -unit_45)),

            # reverse check with std_axes
            (c1, r1, (half_width_r1 + r_c1 - d) .* -i_cap, std_axes, PE2D.Manifold(d, PE2D.rotate_90(std_axes), (r_c1 - d / 2) .* -i_cap)),
            (c1, r1, (half_width_r1 + d) .* -i_cap, std_axes, PE2D.Manifold(r_c1 - d, PE2D.rotate_90(std_axes), (r_c1 - (r_c1 - d) / 2) .* -i_cap)),
            (c1, r1, (half_width_r1 - d) .* -i_cap, std_axes, PE2D.Manifold(r_c1 + d, PE2D.rotate_90(std_axes), (r_c1 - (r_c1 + d) / 2) .* -i_cap)),
            (c1, r1, (half_width_r1 - d) .* i_cap, std_axes, PE2D.Manifold(r_c1 + d, PE2D.rotate_minus_90(std_axes), (r_c1 - (r_c1 + d) / 2) .* i_cap)),
            (c1, r1, (half_width_r1 + d) .* i_cap, std_axes, PE2D.Manifold(r_c1 - d, PE2D.rotate_minus_90(std_axes), (r_c1 - (r_c1 - d) / 2) .* i_cap)),
            (c1, r1, (half_width_r1 + r_c1 - d) .* i_cap, std_axes, PE2D.Manifold(d, PE2D.rotate_minus_90(std_axes), (r_c1 - d / 2) .* i_cap)),

            (c1, r1, (half_height_r1 + r_c1 - d) .* -j_cap, std_axes, PE2D.Manifold(d, PE2D.rotate_180(std_axes), (r_c1 - d / 2) .* -j_cap)),
            (c1, r1, (half_height_r1 + d) .* -j_cap, std_axes, PE2D.Manifold(r_c1 - d, PE2D.rotate_180(std_axes), (r_c1 - (r_c1 - d) / 2) .* -j_cap)),
            (c1, r1, (half_height_r1 - d) .* -j_cap, std_axes, PE2D.Manifold(r_c1 + d, PE2D.rotate_180(std_axes), (r_c1 - (r_c1 + d) / 2) .* -j_cap)),
            (c1, r1, (half_height_r1 - d) .* j_cap, std_axes, PE2D.Manifold(r_c1 + d, std_axes, (r_c1 - (r_c1 + d) / 2) .* j_cap)),
            (c1, r1, (half_height_r1 + d) .* j_cap, std_axes, PE2D.Manifold(r_c1 - d, std_axes, (r_c1 - (r_c1 - d) / 2) .* j_cap)),
            (c1, r1, (half_height_r1 + r_c1 - d) .* j_cap, std_axes, PE2D.Manifold(d, std_axes, (r_c1 - d / 2) .* j_cap)),

            (c1, r1, top_right_r1 .+ (r_c1 - d) .* unit_45, std_axes, PE2D.Manifold(d, PE2D.Axes(-theta_45), (r_c1 - d / 2) .* unit_45)),

            # reverse check with rotated_axes
            (c1, r1, (r_c1 - d) .* unit_45 .+ PE2D.rotate(top_right_r1, rotated_axes), rotated_axes, PE2D.Manifold(d, PE2D.Axes(-theta_45), (r_c1 - d / 2) .* unit_45)),
            ]

            test_manifold_list(manifold_list)
        end

        # @testset "Rect2D vs. Rect2D" begin
            # manifold_list = [
            # # rotated_axes
            # # (r1, r2, (half_height_r1 - d) .* j_cap .+ PE2D.rotate(top_right_r2, rotated_axes), rotated_axes, PE2D.Manifold(d, std_axes, (half_height_r1 - d/2) .* j_cap)),

            # ]

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
