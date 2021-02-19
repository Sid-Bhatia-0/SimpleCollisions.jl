import PhysicsEngine2D
import PhysicsEngine2D: PE2D
import GeometryBasics
const GB = GeometryBasics
import LinearAlgebra
const LA = LinearAlgebra
using Test

function test_collision_list_no_axes(collision_list_no_axes)
    for (a, b, pos_ba, value) in collision_list_no_axes
        @test PE2D.is_colliding(a, b, pos_ba) == value
    end
end

function test_collision_list(collision_list)
    for (a, b, pos_ba, axes_ba, value) in collision_list
        @test PE2D.is_colliding(a, b, pos_ba, axes_ba) == value
    end
end

function test_manifold_list_no_axes(manifold_list_no_axes)
    for (i, (a, b, pos_ba, value)) in enumerate(manifold_list_no_axes)
        manifold_ba = PE2D.Manifold(a, b, pos_ba)
        @show a
        @show b
        @show pos_ba
        @show value
        @show manifold_ba

        @test PE2D.get_penetration(manifold_ba) ≈ PE2D.get_penetration(value)
        @test PE2D.get_normal(manifold_ba) ≈ PE2D.get_normal(value)
        @test PE2D.get_tangent(manifold_ba) ≈ PE2D.get_tangent(value)
        @test PE2D.get_contact(manifold_ba) ≈ PE2D.get_contact(value)
    end
end

function test_manifold_list(manifold_list)
    for (i, (a, b, pos_ba, axes_ba, value)) in enumerate(manifold_list)
        manifold_ba = PE2D.Manifold(a, b, pos_ba, axes_ba)
        @test PE2D.get_penetration(manifold_ba) ≈ PE2D.get_penetration(value)
        @test PE2D.get_normal(manifold_ba) ≈ PE2D.get_normal(value)
        @test PE2D.get_tangent(manifold_ba) ≈ PE2D.get_tangent(value)
        @test PE2D.get_contact(manifold_ba) ≈ PE2D.get_contact(value)
    end
end

@testset "PhysicsEngine2D.jl" begin
    T = Float32
    VecType = GB.Vec{2, T}

    origin = zero(VecType)
    d = convert(T, 0.01)

    std_axes = PE2D.Axes{T}()
    i_cap = PE2D.get_x_cap(std_axes)
    j_cap = PE2D.get_y_cap(std_axes)

    theta = convert(T, π / 6)
    rotated_axes = PE2D.Axes(theta)
    theta_45 = convert(T, π / 4)
    unit_45 = (i_cap .+ j_cap) / convert(T, sqrt(2))

    point = PE2D.StdPoint{T}()

    half_length_l1 = one(T)
    l1 = PE2D.StdLine(half_length_l1)
    p1_l1 = PE2D.get_tail(l1)
    p2_l1 = PE2D.get_head(l1)

    half_length_l2 = convert(T, 2)
    l2 = PE2D.StdLine(half_length_l2)
    p1_l2 = PE2D.get_tail(l2)
    p2_l2 = PE2D.get_head(l2)

    r_c1 = one(T)
    c1 = PE2D.StdCircle(r_c1)

    r_c2 = convert(T, 2)
    c2 = PE2D.StdCircle(r_c2)

    half_width_r1 = one(T)
    half_height_r1 = convert(T, 0.5)
    r1 = PE2D.StdRect(half_width_r1, half_height_r1)
    top_right_r1 = PE2D.get_top_right(r1)
    theta_r1 = atan(half_height_r1, half_width_r1)

    half_width_r2 = convert(T, 2)
    half_height_r2 = one(T)
    r2 = PE2D.StdRect(half_width_r2, half_height_r2)
    top_right_r2 = PE2D.get_top_right(r2)
    theta_r2 = atan(half_height_r2, half_width_r2)

    @testset "Area" begin
        @testset "StdRect" begin
            @test PE2D.get_area(r2) == convert(T, 8)
        end

        @testset "StdCircle" begin
            @test PE2D.get_area(c1) ≈ convert(T, π)
        end
    end

    @testset "Collision detection" begin
        @testset "StdLine vs. StdLine" begin
            collision_list = [
            # std_axes
            (l1, l2, origin, std_axes, true),

            (l1, l2, (half_length_l1 + half_length_l2 + d) .* -i_cap, std_axes, false),
            (l1, l2, (half_length_l1 + half_length_l2 - d) .* -i_cap, std_axes, true),
            (l1, l2, (half_length_l1 + half_length_l2 - d) .* i_cap, std_axes, true),
            (l1, l2, (half_length_l1 + half_length_l2 + d) .* i_cap, std_axes, false),

            (l1, l2, d .* -j_cap, std_axes, false),
            (l1, l2, d .* j_cap, std_axes, false),

            # rotated_axes
            (l1, l2, origin, rotated_axes, true),

            (l1, l2, (half_length_l1 + half_length_l2 + d) .* -i_cap, rotated_axes, false),
            (l1, l2, (half_length_l1 + half_length_l2 - d) .* -i_cap, rotated_axes, false),
            (l1, l2, (half_length_l1 + d) .* -i_cap, rotated_axes, false),
            (l1, l2, (half_length_l1 - d) .* -i_cap, rotated_axes, true),
            (l1, l2, (half_length_l1 - d) .* i_cap, rotated_axes, true),
            (l1, l2, (half_length_l1 + d) .* i_cap, rotated_axes, false),
            (l1, l2, (half_length_l1 + half_length_l2 - d) .* i_cap, rotated_axes, false),
            (l1, l2, (half_length_l1 + half_length_l2 + d) .* i_cap, rotated_axes, false),

            (l2, l1, (half_length_l1 * sin(theta) + d) .* -j_cap, rotated_axes, false),
            (l2, l1, (half_length_l1 * sin(theta) - d) .* -j_cap, rotated_axes, true),
            (l2, l1, d .* -j_cap, rotated_axes, true),
            (l2, l1, d .* j_cap, rotated_axes, true),
            (l2, l1, (half_length_l1 * sin(theta) - d) .* j_cap, rotated_axes, true),
            (l2, l1, (half_length_l1 * sin(theta) + d) .* j_cap, rotated_axes, false),
            ]

            test_collision_list(collision_list)

            collision_list_no_axes = [
            # std_axes
            (l1, l2, origin, true),

            (l1, l2, (half_length_l1 + half_length_l2 + d) .* -i_cap, false),
            (l1, l2, (half_length_l1 + half_length_l2 - d) .* -i_cap, true),
            (l1, l2, (half_length_l1 + half_length_l2 - d) .* i_cap, true),
            (l1, l2, (half_length_l1 + half_length_l2 + d) .* i_cap, false),

            (l1, l2, d .* -j_cap, false),
            (l1, l2, d .* j_cap, false),
            ]

            test_collision_list_no_axes(collision_list_no_axes)
        end

        @testset "StdCircle vs. StdPoint" begin
            collision_list = [
            # std_axes
            (c1, point, origin, std_axes, true),

            (c1, point, (r_c1 + d) .* -i_cap, std_axes, false),
            (c1, point, (r_c1 - d) .* -i_cap, std_axes, true),
            (c1, point, (r_c1 - d) .* i_cap, std_axes, true),
            (c1, point, (r_c1 + d) .* i_cap, std_axes, false),

            (c1, point, (r_c1 + d) .* -j_cap, std_axes, false),
            (c1, point, (r_c1 - d) .* -j_cap, std_axes, true),
            (c1, point, (r_c1 - d) .* j_cap, std_axes, true),
            (c1, point, (r_c1 + d) .* j_cap, std_axes, false),

            (c1, point, (r_c1 + d) .* unit_45, std_axes, false),
            (c1, point, (r_c1 - d) .* unit_45, std_axes, true),

            # reverse check with std_axes
            (point, c1, origin, std_axes, true),

            (point, c1, (r_c1 + d) .* -i_cap, std_axes, false),
            (point, c1, (r_c1 - d) .* -i_cap, std_axes, true),
            (point, c1, (r_c1 - d) .* i_cap, std_axes, true),
            (point, c1, (r_c1 + d) .* i_cap, std_axes, false),

            (point, c1, (r_c1 + d) .* -j_cap, std_axes, false),
            (point, c1, (r_c1 - d) .* -j_cap, std_axes, true),
            (point, c1, (r_c1 - d) .* j_cap, std_axes, true),
            (point, c1, (r_c1 + d) .* j_cap, std_axes, false),

            (point, c1, (r_c1 + d) .* unit_45, std_axes, false),
            (point, c1, (r_c1 - d) .* unit_45, std_axes, true),
            ]

            test_collision_list(collision_list)

            collision_list_no_axes = [
            # std_axes
            (c1, point, origin, true),

            (c1, point, (r_c1 + d) .* -i_cap, false),
            (c1, point, (r_c1 - d) .* -i_cap, true),
            (c1, point, (r_c1 - d) .* i_cap, true),
            (c1, point, (r_c1 + d) .* i_cap, false),

            (c1, point, (r_c1 + d) .* -j_cap, false),
            (c1, point, (r_c1 - d) .* -j_cap, true),
            (c1, point, (r_c1 - d) .* j_cap, true),
            (c1, point, (r_c1 + d) .* j_cap, false),

            (c1, point, (r_c1 + d) .* unit_45, false),
            (c1, point, (r_c1 - d) .* unit_45, true),

            # reverse check with std_axes
            (point, c1, origin, true),

            (point, c1, (r_c1 + d) .* -i_cap, false),
            (point, c1, (r_c1 - d) .* -i_cap, true),
            (point, c1, (r_c1 - d) .* i_cap, true),
            (point, c1, (r_c1 + d) .* i_cap, false),

            (point, c1, (r_c1 + d) .* -j_cap, false),
            (point, c1, (r_c1 - d) .* -j_cap, true),
            (point, c1, (r_c1 - d) .* j_cap, true),
            (point, c1, (r_c1 + d) .* j_cap, false),

            (point, c1, (r_c1 + d) .* unit_45, false),
            (point, c1, (r_c1 - d) .* unit_45, true),
            ]

            test_collision_list_no_axes(collision_list_no_axes)
        end

        @testset "StdCircle vs. StdLine" begin
            collision_list = [
            # std_axes
            (l1, c2, origin, std_axes, true),
            (l2, c1, origin, std_axes, true),

            (l2, c1, (half_length_l2 + r_c1 + d) .* -i_cap, std_axes, false),
            (l2, c1, (half_length_l2 + r_c1 - d) .* -i_cap, std_axes, true),
            (l2, c1, (half_length_l2 + r_c1 - d) .* i_cap, std_axes, true),
            (l2, c1, (half_length_l2 + r_c1 + d) .* i_cap, std_axes, false),

            (l2, c1, (r_c1 + d) .* -j_cap, std_axes, false),
            (l2, c1, (r_c1 - d) .* -j_cap, std_axes, true),
            (l2, c1, (r_c1 - d) .* j_cap, std_axes, true),
            (l2, c1, (r_c1 + d) .* j_cap, std_axes, false),

            # reverse check with std_axes
            (c2, l1, origin, std_axes, true),
            (c1, l2, origin, std_axes, true),

            (c1, l2, (half_length_l2 + r_c1 + d) .* -i_cap, std_axes, false),
            (c1, l2, (half_length_l2 + r_c1 - d) .* -i_cap, std_axes, true),
            (c1, l2, (half_length_l2 + r_c1 - d) .* i_cap, std_axes, true),
            (c1, l2, (half_length_l2 + r_c1 + d) .* i_cap, std_axes, false),

            (c1, l2, (r_c1 + d) .* -j_cap, std_axes, false),
            (c1, l2, (r_c1 - d) .* -j_cap, std_axes, true),
            (c1, l2, (r_c1 - d) .* j_cap, std_axes, true),
            (c1, l2, (r_c1 + d) .* j_cap, std_axes, false),

            # reverse check with rotated_axes
            (c2, l1, origin, rotated_axes, true),
            (c1, l2, origin, rotated_axes, true),

            (c1, l1, (sqrt(r_c1 ^ 2 - (half_length_l1 * sin(theta)) ^ 2) + half_length_l1 * cos(theta) + d) .* -i_cap, rotated_axes, false),
            (c1, l1, (sqrt(r_c1 ^ 2 - (half_length_l1 * sin(theta)) ^ 2) + half_length_l1 * cos(theta) - d) .* -i_cap, rotated_axes, true),
            (c1, l1, (sqrt(r_c1 ^ 2 - (half_length_l1 * sin(theta)) ^ 2) + half_length_l1 * cos(theta) - d) .* i_cap, rotated_axes, true),
            (c1, l1, (sqrt(r_c1 ^ 2 - (half_length_l1 * sin(theta)) ^ 2) + half_length_l1 * cos(theta) + d) .* i_cap, rotated_axes, false),

            (c1, l2, (r_c1 / cos(theta) + d) .* -j_cap, rotated_axes, false),
            (c1, l2, (r_c1 / cos(theta) - d) .* -j_cap, rotated_axes, true),
            (c1, l2, (r_c1 / cos(theta) - d) .* j_cap, rotated_axes, true),
            (c1, l2, (r_c1 / cos(theta) + d) .* j_cap, rotated_axes, false),
            ]

            test_collision_list(collision_list)

            collision_list_no_axes = [
            # std_axes
            (l1, c2, origin, true),
            (l2, c1, origin, true),

            (l2, c1, (half_length_l2 + r_c1 + d) .* -i_cap, false),
            (l2, c1, (half_length_l2 + r_c1 - d) .* -i_cap, true),
            (l2, c1, (half_length_l2 + r_c1 - d) .* i_cap, true),
            (l2, c1, (half_length_l2 + r_c1 + d) .* i_cap, false),

            (l2, c1, (r_c1 + d) .* -j_cap, false),
            (l2, c1, (r_c1 - d) .* -j_cap, true),
            (l2, c1, (r_c1 - d) .* j_cap, true),
            (l2, c1, (r_c1 + d) .* j_cap, false),

            # reverse check with std_axes
            (c2, l1, origin, true),
            (c1, l2, origin, true),

            (c1, l2, (half_length_l2 + r_c1 + d) .* -i_cap, false),
            (c1, l2, (half_length_l2 + r_c1 - d) .* -i_cap, true),
            (c1, l2, (half_length_l2 + r_c1 - d) .* i_cap, true),
            (c1, l2, (half_length_l2 + r_c1 + d) .* i_cap, false),

            (c1, l2, (r_c1 + d) .* -j_cap, false),
            (c1, l2, (r_c1 - d) .* -j_cap, true),
            (c1, l2, (r_c1 - d) .* j_cap, true),
            (c1, l2, (r_c1 + d) .* j_cap, false),
            ]

            test_collision_list_no_axes(collision_list_no_axes)
        end

        @testset "StdCircle vs. StdCircle" begin
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

            collision_list_no_axes = [
            # std_axes
            (c1, c2, origin, true),

            (c1, c2, (r_c1 + r_c2 + d) .* -i_cap, false),
            (c1, c2, (r_c1 + r_c2 - d) .* -i_cap, true),
            (c1, c2, (r_c1 + r_c2 - d) .* i_cap, true),
            (c1, c2, (r_c1 + r_c2 + d) .* i_cap, false),

            (c1, c2, (r_c1 + r_c2 + d) .* -j_cap, false),
            (c1, c2, (r_c1 + r_c2 - d) .* -j_cap, true),
            (c1, c2, (r_c1 + r_c2 - d) .* j_cap, true),
            (c1, c2, (r_c1 + r_c2 + d) .* j_cap, false),

            (c1, c2, (r_c1 + r_c2 + d) .* unit_45, false),
            (c1, c2, (r_c1 + r_c2 - d) .* unit_45, true),
            ]

            test_collision_list_no_axes(collision_list_no_axes)
        end

        @testset "StdRect vs. StdPoint" begin
            collision_list = [
            # std_axes
            (r1, point, origin, std_axes, true),

            (r1, point, (half_width_r1 + d) .* -i_cap, std_axes, false),
            (r1, point, (half_width_r1 - d) .* -i_cap, std_axes, true),
            (r1, point, (half_width_r1 - d) .* i_cap, std_axes, true),
            (r1, point, (half_width_r1 + d) .* i_cap, std_axes, false),

            (r1, point, (half_height_r1 + d) .* -j_cap, std_axes, false),
            (r1, point, (half_height_r1 - d) .* -j_cap, std_axes, true),
            (r1, point, (half_height_r1 - d) .* j_cap, std_axes, true),
            (r1, point, (half_height_r1 + d) .* j_cap, std_axes, false),

            (r1, point, top_right_r1 .+ d, std_axes, false),
            (r1, point, top_right_r1 .- d, std_axes, true),

            # reverse check with std_axes
            (point, r1, origin, std_axes, true),

            (point, r1, (half_width_r1 + d) .* -i_cap, std_axes, false),
            (point, r1, (half_width_r1 - d) .* -i_cap, std_axes, true),
            (point, r1, (half_width_r1 - d) .* i_cap, std_axes, true),
            (point, r1, (half_width_r1 + d) .* i_cap, std_axes, false),

            (point, r1, (half_height_r1 + d) .* -j_cap, std_axes, false),
            (point, r1, (half_height_r1 - d) .* -j_cap, std_axes, true),
            (point, r1, (half_height_r1 - d) .* j_cap, std_axes, true),
            (point, r1, (half_height_r1 + d) .* j_cap, std_axes, false),

            (point, r1, top_right_r1 .+ d, std_axes, false),
            (point, r1, top_right_r1 .- d, std_axes, true),

            # reverse check with rotated_axes
            (point, r1, origin, rotated_axes, true),

            (point, r1, (half_height_r1 / sin(theta) + d) .* -i_cap, rotated_axes, false),
            (point, r1, (half_height_r1 / sin(theta) - d) .* -i_cap, rotated_axes, true),
            (point, r1, (half_height_r1 / sin(theta) - d) .* i_cap, rotated_axes, true),
            (point, r1, (half_height_r1 / sin(theta) + d) .* i_cap, rotated_axes, false),

            (point, r1, (half_height_r1 / cos(theta) + d) .* -j_cap, rotated_axes, false),
            (point, r1, (half_height_r1 / cos(theta) - d) .* -j_cap, rotated_axes, true),
            (point, r1, (half_height_r1 / cos(theta) - d) .* j_cap, rotated_axes, true),
            (point, r1, (half_height_r1 / cos(theta) + d) .* j_cap, rotated_axes, false),
            ]

            test_collision_list(collision_list)

            collision_list_no_axes = [
            # std_axes
            (r1, point, origin, true),

            (r1, point, (half_width_r1 + d) .* -i_cap, false),
            (r1, point, (half_width_r1 - d) .* -i_cap, true),
            (r1, point, (half_width_r1 - d) .* i_cap, true),
            (r1, point, (half_width_r1 + d) .* i_cap, false),

            (r1, point, (half_height_r1 + d) .* -j_cap, false),
            (r1, point, (half_height_r1 - d) .* -j_cap, true),
            (r1, point, (half_height_r1 - d) .* j_cap, true),
            (r1, point, (half_height_r1 + d) .* j_cap, false),

            (r1, point, top_right_r1 .+ d, false),
            (r1, point, top_right_r1 .- d, true),

            # reverse check with std_axes
            (point, r1, origin, true),

            (point, r1, (half_width_r1 + d) .* -i_cap, false),
            (point, r1, (half_width_r1 - d) .* -i_cap, true),
            (point, r1, (half_width_r1 - d) .* i_cap, true),
            (point, r1, (half_width_r1 + d) .* i_cap, false),

            (point, r1, (half_height_r1 + d) .* -j_cap, false),
            (point, r1, (half_height_r1 - d) .* -j_cap, true),
            (point, r1, (half_height_r1 - d) .* j_cap, true),
            (point, r1, (half_height_r1 + d) .* j_cap, false),

            (point, r1, top_right_r1 .+ d, false),
            (point, r1, top_right_r1 .- d, true),
            ]

            test_collision_list_no_axes(collision_list_no_axes)
        end

        @testset "StdRect vs. StdLine" begin
            collision_list = [
            # std_axes
            (r1, l1, origin, std_axes, true),
            (r1, l2, origin, std_axes, true),

            (r1, l1, (half_width_r1 + half_length_l1 + d) .* -i_cap, std_axes, false),
            (r1, l1, (half_width_r1 + half_length_l1 - d) .* -i_cap, std_axes, true),
            (r1, l1, (half_width_r1 + half_length_l1 - d) .* i_cap, std_axes, true),
            (r1, l1, (half_width_r1 + half_length_l1 + d) .* i_cap, std_axes, false),

            (r1, l1, (half_height_r1 + d) .* -j_cap, std_axes, false),
            (r1, l1, (half_height_r1 - d) .* -j_cap, std_axes, true),
            (r1, l1, (half_height_r1 - d) .* j_cap, std_axes, true),
            (r1, l1, (half_height_r1 + d) .* j_cap, std_axes, false),

            # rotated_axes
            (r1, l1, origin, rotated_axes, true),
            (r1, l2, origin, rotated_axes, true),

            (r2, l1, (half_width_r2 + half_length_l1 * cos(theta) + d) .* -i_cap, rotated_axes, false),
            (r2, l1, (half_width_r2 + half_length_l1 * cos(theta) - d) .* -i_cap, rotated_axes, true),
            (r2, l1, (half_width_r2 + half_length_l1 * cos(theta) - d) .* i_cap, rotated_axes, true),
            (r2, l1, (half_width_r2 + half_length_l1 * cos(theta) + d) .* i_cap, rotated_axes, false),

            (r2, l1, (half_height_r2 + half_length_l1 * sin(theta) + d) .* -j_cap, rotated_axes, false),
            (r2, l1, (half_height_r2 + half_length_l1 * sin(theta) - d) .* -j_cap, rotated_axes, true),
            (r2, l1, (half_height_r2 + half_length_l1 * sin(theta) - d) .* j_cap, rotated_axes, true),
            (r2, l1, (half_height_r2 + half_length_l1 * sin(theta) + d) .* j_cap, rotated_axes, false),

            # reverse check with std_axes
            (l1, r1, origin, std_axes, true),
            (l2, r1, origin, std_axes, true),

            (l1, r1, (half_width_r1 + half_length_l1 + d) .* -i_cap, std_axes, false),
            (l1, r1, (half_width_r1 + half_length_l1 - d) .* -i_cap, std_axes, true),
            (l1, r1, (half_width_r1 + half_length_l1 - d) .* i_cap, std_axes, true),
            (l1, r1, (half_width_r1 + half_length_l1 + d) .* i_cap, std_axes, false),

            (l1, r1, (half_height_r1 + d) .* -j_cap, std_axes, false),
            (l1, r1, (half_height_r1 - d) .* -j_cap, std_axes, true),
            (l1, r1, (half_height_r1 - d) .* j_cap, std_axes, true),
            (l1, r1, (half_height_r1 + d) .* j_cap, std_axes, false),

            # reverse check with rotated_axes
            (l1, r1, origin, rotated_axes, true),
            (l2, r1, origin, rotated_axes, true),

            (l2, r2, (half_height_r2 / sin(theta) + half_length_l2 + d) .* -i_cap, rotated_axes, false),
            (l2, r2, (half_height_r2 / sin(theta) + half_length_l2 - d) .* -i_cap, rotated_axes, true),
            (l2, r2, (half_height_r2 / sin(theta) + half_length_l2 - d) .* i_cap, rotated_axes, true),
            (l2, r2, (half_height_r2 / sin(theta) + half_length_l2 + d) .* i_cap, rotated_axes, false),

            (l2, r2, (half_width_r2 * sin(theta) + half_height_r2 * cos(theta) + d) .* -j_cap, rotated_axes, false),
            (l2, r2, (half_width_r2 * sin(theta) + half_height_r2 * cos(theta) - d) .* -j_cap, rotated_axes, true),
            (l2, r2, (half_width_r2 * sin(theta) + half_height_r2 * cos(theta) - d) .* j_cap, rotated_axes, true),
            (l2, r2, (half_width_r2 * sin(theta) + half_height_r2 * cos(theta) + d) .* j_cap, rotated_axes, false),
            ]

            test_collision_list(collision_list)

            collision_list_no_axes = [
            # std_axes
            (r1, l1, origin, true),
            (r1, l2, origin, true),

            (r1, l1, (half_width_r1 + half_length_l1 + d) .* -i_cap, false),
            (r1, l1, (half_width_r1 + half_length_l1 - d) .* -i_cap, true),
            (r1, l1, (half_width_r1 + half_length_l1 - d) .* i_cap, true),
            (r1, l1, (half_width_r1 + half_length_l1 + d) .* i_cap, false),

            (r1, l1, (half_height_r1 + d) .* -j_cap, false),
            (r1, l1, (half_height_r1 - d) .* -j_cap, true),
            (r1, l1, (half_height_r1 - d) .* j_cap, true),
            (r1, l1, (half_height_r1 + d) .* j_cap, false),

            # reverse check with std_axes
            (l1, r1, origin, true),
            (l2, r1, origin, true),

            (l1, r1, (half_width_r1 + half_length_l1 + d) .* -i_cap, false),
            (l1, r1, (half_width_r1 + half_length_l1 - d) .* -i_cap, true),
            (l1, r1, (half_width_r1 + half_length_l1 - d) .* i_cap, true),
            (l1, r1, (half_width_r1 + half_length_l1 + d) .* i_cap, false),

            (l1, r1, (half_height_r1 + d) .* -j_cap, false),
            (l1, r1, (half_height_r1 - d) .* -j_cap, true),
            (l1, r1, (half_height_r1 - d) .* j_cap, true),
            (l1, r1, (half_height_r1 + d) .* j_cap, false),
            ]

            test_collision_list_no_axes(collision_list_no_axes)
        end

        @testset "StdRect vs. StdCircle" begin
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

            collision_list_no_axes = [
            # std_axes
            (r1, c1, origin, true),

            (r1, c1, (half_width_r1 + r_c1 + d) .* -i_cap, false),
            (r1, c1, (half_width_r1 + r_c1 - d) .* -i_cap, true),
            (r1, c1, (half_width_r1 + r_c1 - d) .* i_cap, true),
            (r1, c1, (half_width_r1 + r_c1 + d) .* i_cap, false),

            (r1, c1, (half_height_r1 + r_c1 + d) .* -j_cap, false),
            (r1, c1, (half_height_r1 + r_c1 - d) .* -j_cap, true),
            (r1, c1, (half_height_r1 + r_c1 - d) .* j_cap, true),
            (r1, c1, (half_height_r1 + r_c1 + d) .* -j_cap, false),

            (r1, c1, top_right_r1 .+ (r_c1 + d) .* unit_45, false),
            (r1, c1, top_right_r1 .+ (r_c1 - d) .* unit_45, true),

            # reverse check with std_axes
            (c1, r1, origin, true),

            (c1, r1, (half_width_r1 + r_c1 + d) .* -i_cap, false),
            (c1, r1, (half_width_r1 + r_c1 - d) .* -i_cap, true),
            (c1, r1, (half_width_r1 + r_c1 - d) .* i_cap, true),
            (c1, r1, (half_width_r1 + r_c1 + d) .* i_cap, false),

            (c1, r1, (half_height_r1 + r_c1 + d) .* -j_cap, false),
            (c1, r1, (half_height_r1 + r_c1 - d) .* -j_cap, true),
            (c1, r1, (half_height_r1 + r_c1 - d) .* j_cap, true),
            (c1, r1, (half_height_r1 + r_c1 + d) .* -j_cap, false),

            (c1, r1, top_right_r1 .+ (r_c1 + d) .* unit_45, false),
            (c1, r1, top_right_r1 .+ (r_c1 - d) .* unit_45, true),
            ]

            test_collision_list_no_axes(collision_list_no_axes)
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

            collision_list_no_axes = [
            # std_axes
            (r2, r1, origin, true),

            (r2, r1, (half_width_r1 + half_width_r2 + d) .* -i_cap, false),
            (r2, r1, (half_width_r1 + half_width_r2 - d) .* -i_cap, true),
            (r2, r1, (half_width_r1 + half_width_r2 - d) .* i_cap, true),
            (r2, r1, (half_width_r1 + half_width_r2 + d) .* i_cap, false),

            (r2, r1, (half_height_r1 + half_height_r2 + d) .* -j_cap, false),
            (r2, r1, (half_height_r1 + half_height_r2 - d) .* -j_cap, true),
            (r2, r1, (half_height_r1 + half_height_r2 - d) .* j_cap, true),
            (r2, r1, (half_height_r1 + half_height_r2 + d) .* j_cap, false),

            (r2, r1, top_right_r1 .+ top_right_r2 .+ d, false),
            (r2, r1, top_right_r1 .+ top_right_r2 .- d, true),
            ]

            test_collision_list_no_axes(collision_list_no_axes)
        end
    end

    @testset "Manifold generation" begin
        @testset "StdCircle vs. StdCircle" begin
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

        @testset "StdRect vs. StdCircle" begin
            manifold_list_no_axes = [
            # std_axes
            (r1, c1, (half_width_r1 + r_c1 - d) .* -i_cap, PE2D.Manifold(d, PE2D.rotate_90(std_axes), (half_width_r1 - d / 2) .* -i_cap)),
            (r1, c1, (half_width_r1 + d) .* -i_cap, PE2D.Manifold(r_c1 - d, PE2D.rotate_90(std_axes), (half_width_r1 - (r_c1 - d) / 2) .* -i_cap)),
            (r1, c1, (half_width_r1 - d) .* -i_cap, PE2D.Manifold(r_c1 + d, PE2D.rotate_90(std_axes), (half_width_r1 - (r_c1 + d) / 2) .* -i_cap)),
            (r1, c1, (half_width_r1 - d) .* i_cap, PE2D.Manifold(r_c1 + d, PE2D.rotate_minus_90(std_axes), (half_width_r1 - (r_c1 + d) / 2) .* i_cap)),
            (r1, c1, (half_width_r1 + d) .* i_cap, PE2D.Manifold(r_c1 - d, PE2D.rotate_minus_90(std_axes), (half_width_r1 - (r_c1 - d) / 2) .* i_cap)),
            (r1, c1, (half_width_r1 + r_c1 - d) .* i_cap, PE2D.Manifold(d, PE2D.rotate_minus_90(std_axes), (half_width_r1 - d / 2) .* i_cap)),

            (r1, c1, (half_height_r1 + r_c1 - d) .* -j_cap, PE2D.Manifold(d, PE2D.rotate_180(std_axes), (half_height_r1 - d / 2) .* -j_cap)),
            (r1, c1, (half_height_r1 + d) .* -j_cap, PE2D.Manifold(r_c1 - d, PE2D.rotate_180(std_axes), (half_height_r1 - (r_c1 - d) / 2) .* -j_cap)),
            (r1, c1, (half_height_r1 - d) .* -j_cap, PE2D.Manifold(r_c1 + d, PE2D.rotate_180(std_axes), (half_height_r1 - (r_c1 + d) / 2) .* -j_cap)),
            (r1, c1, (half_height_r1 - d) .* j_cap, PE2D.Manifold(r_c1 + d, std_axes, (half_height_r1 - (r_c1 + d) / 2) .* j_cap)),
            (r1, c1, (half_height_r1 + d) .* j_cap, PE2D.Manifold(r_c1 - d, std_axes, (half_height_r1 - (r_c1 - d) / 2) .* j_cap)),
            (r1, c1, (half_height_r1 + r_c1 - d) .* j_cap, PE2D.Manifold(d, std_axes, (half_height_r1 - d / 2) .* j_cap)),

            (r1, c1, top_right_r1 .+ (r_c1 - d) .* unit_45, PE2D.Manifold(d, PE2D.Axes(-theta_45), top_right_r1 .+ (d / 2) .* -unit_45)),

            # reverse check with std_axes
            (c1, r1, (half_width_r1 + r_c1 - d) .* -i_cap, PE2D.Manifold(d, PE2D.rotate_90(std_axes), (r_c1 - d / 2) .* -i_cap)),
            (c1, r1, (half_width_r1 + d) .* -i_cap, PE2D.Manifold(r_c1 - d, PE2D.rotate_90(std_axes), (r_c1 - (r_c1 - d) / 2) .* -i_cap)),
            (c1, r1, (half_width_r1 - d) .* -i_cap, PE2D.Manifold(r_c1 + d, PE2D.rotate_90(std_axes), (r_c1 - (r_c1 + d) / 2) .* -i_cap)),
            (c1, r1, (half_width_r1 - d) .* i_cap, PE2D.Manifold(r_c1 + d, PE2D.rotate_minus_90(std_axes), (r_c1 - (r_c1 + d) / 2) .* i_cap)),
            (c1, r1, (half_width_r1 + d) .* i_cap, PE2D.Manifold(r_c1 - d, PE2D.rotate_minus_90(std_axes), (r_c1 - (r_c1 - d) / 2) .* i_cap)),
            (c1, r1, (half_width_r1 + r_c1 - d) .* i_cap, PE2D.Manifold(d, PE2D.rotate_minus_90(std_axes), (r_c1 - d / 2) .* i_cap)),

            (c1, r1, (half_height_r1 + r_c1 - d) .* -j_cap, PE2D.Manifold(d, PE2D.rotate_180(std_axes), (r_c1 - d / 2) .* -j_cap)),
            (c1, r1, (half_height_r1 + d) .* -j_cap, PE2D.Manifold(r_c1 - d, PE2D.rotate_180(std_axes), (r_c1 - (r_c1 - d) / 2) .* -j_cap)),
            (c1, r1, (half_height_r1 - d) .* -j_cap, PE2D.Manifold(r_c1 + d, PE2D.rotate_180(std_axes), (r_c1 - (r_c1 + d) / 2) .* -j_cap)),
            (c1, r1, (half_height_r1 - d) .* j_cap, PE2D.Manifold(r_c1 + d, std_axes, (r_c1 - (r_c1 + d) / 2) .* j_cap)),
            (c1, r1, (half_height_r1 + d) .* j_cap, PE2D.Manifold(r_c1 - d, std_axes, (r_c1 - (r_c1 - d) / 2) .* j_cap)),
            (c1, r1, (half_height_r1 + r_c1 - d) .* j_cap, PE2D.Manifold(d, std_axes, (r_c1 - d / 2) .* j_cap)),

            (c1, r1, top_right_r1 .+ (r_c1 - d) .* unit_45, PE2D.Manifold(d, PE2D.Axes(-theta_45), (r_c1 - d / 2) .* unit_45)),
            ]

            test_manifold_list_no_axes(manifold_list_no_axes)

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

        @testset "Rect2D vs. Rect2D" begin
            manifold_list_no_axes = [
            # std_axes
            (r2, r1, (half_height_r2 + half_height_r1 - d) .* -j_cap, PE2D.Manifold(d, PE2D.rotate_180(std_axes), (half_height_r2 - d/2) .* -j_cap)),
            (r2, r1, (half_height_r2 + d) .* -j_cap, PE2D.Manifold(half_height_r1 - d, PE2D.rotate_180(std_axes), (half_height_r2 - (half_height_r1 - d)/2) .* -j_cap)),
            (r2, r1, (half_height_r2 - d) .* -j_cap, PE2D.Manifold(half_height_r1 + d, PE2D.rotate_180(std_axes), (half_height_r2 - (half_height_r1 + d)/2) .* -j_cap)),
            (r2, r1, d .* -j_cap, PE2D.Manifold(half_height_r2 + half_height_r1 - d, PE2D.rotate_180(std_axes), d .* -j_cap)),
            (r2, r1, d .* j_cap, PE2D.Manifold(half_height_r2 + half_height_r1 - d, std_axes, d .* j_cap)),
            (r2, r1, (half_height_r2 - d) .* j_cap, PE2D.Manifold(half_height_r1 + d, std_axes, (half_height_r2 - (half_height_r1 + d)/2) .* j_cap)),
            (r2, r1, (half_height_r2 + d) .* j_cap, PE2D.Manifold(half_height_r1 - d, std_axes, (half_height_r2 - (half_height_r1 - d)/2) .* j_cap)),
            (r2, r1, (half_height_r2 + half_height_r1 - d) .* j_cap, PE2D.Manifold(d, std_axes, (half_height_r2 - d/2) .* j_cap)),

            (r2, r1, (half_width_r2 + half_width_r1 - d) .* -i_cap, PE2D.Manifold(d, PE2D.rotate_90(std_axes), (half_width_r2 - d/2) .* -i_cap)),
            (r2, r1, (half_width_r2 + d) .* -i_cap, PE2D.Manifold(half_width_r1 - d, PE2D.rotate_90(std_axes), (half_width_r2 - (half_width_r1 - d)/2) .* -i_cap)),
            (r2, r1, (half_width_r2 - d) .* -i_cap, PE2D.Manifold(half_width_r1 + d, PE2D.rotate_90(std_axes), (half_width_r2 - (half_width_r1 + d)/2) .* -i_cap)),
            (r2, r1, (half_width_r2 - d) .* i_cap, PE2D.Manifold(half_width_r1 + d, PE2D.rotate_minus_90(std_axes), (half_width_r2 - (half_width_r1 + d)/2) .* i_cap)),
            (r2, r1, (half_width_r2 + d) .* i_cap, PE2D.Manifold(half_width_r1 - d, PE2D.rotate_minus_90(std_axes), (half_width_r2 - (half_width_r1 - d)/2) .* i_cap)),
            (r2, r1, (half_width_r2 + half_width_r1 - d) .* i_cap, PE2D.Manifold(d, PE2D.rotate_minus_90(std_axes), (half_width_r2 - d/2) .* i_cap)),

            (r2, r1, top_right_r2 .- d, PE2D.Manifold(half_height_r1 + d, std_axes, top_right_r2 .+ (half_width_r1 + d)/2 .* -i_cap .+ (half_height_r1 + d)/2 .* -j_cap)),
            (r2, r1, top_right_r2, PE2D.Manifold(half_height_r1, std_axes, top_right_r2 .+ (half_width_r1/2) .* -i_cap .+ (half_height_r1/2) .* -j_cap)),
            (r2, r1, top_right_r2 .+ d, PE2D.Manifold(half_height_r1 - d, std_axes, top_right_r2 .+ (half_width_r1 - d)/2 .* -i_cap .+ (half_height_r1 - d)/2 .* -j_cap)),
            ]

            test_manifold_list_no_axes(manifold_list_no_axes)

            manifold_list = [
            # std_axes
            (r2, r1, (half_height_r2 + half_height_r1 - d) .* -j_cap, std_axes, PE2D.Manifold(d, PE2D.rotate_180(std_axes), (half_height_r2 - d/2) .* -j_cap)),
            (r2, r1, (half_height_r2 + d) .* -j_cap, std_axes, PE2D.Manifold(half_height_r1 - d, PE2D.rotate_180(std_axes), (half_height_r2 - (half_height_r1 - d)/2) .* -j_cap)),
            (r2, r1, (half_height_r2 - d) .* -j_cap, std_axes, PE2D.Manifold(half_height_r1 + d, PE2D.rotate_180(std_axes), (half_height_r2 - (half_height_r1 + d)/2) .* -j_cap)),
            (r2, r1, d .* -j_cap, std_axes, PE2D.Manifold(half_height_r2 + half_height_r1 - d, PE2D.rotate_180(std_axes), d .* -j_cap)),
            (r2, r1, d .* j_cap, std_axes, PE2D.Manifold(half_height_r2 + half_height_r1 - d, std_axes, d .* j_cap)),
            (r2, r1, (half_height_r2 - d) .* j_cap, std_axes, PE2D.Manifold(half_height_r1 + d, std_axes, (half_height_r2 - (half_height_r1 + d)/2) .* j_cap)),
            (r2, r1, (half_height_r2 + d) .* j_cap, std_axes, PE2D.Manifold(half_height_r1 - d, std_axes, (half_height_r2 - (half_height_r1 - d)/2) .* j_cap)),
            (r2, r1, (half_height_r2 + half_height_r1 - d) .* j_cap, std_axes, PE2D.Manifold(d, std_axes, (half_height_r2 - d/2) .* j_cap)),

            (r2, r1, (half_width_r2 + half_width_r1 - d) .* -i_cap, std_axes, PE2D.Manifold(d, PE2D.rotate_90(std_axes), (half_width_r2 - d/2) .* -i_cap)),
            (r2, r1, (half_width_r2 + d) .* -i_cap, std_axes, PE2D.Manifold(half_width_r1 - d, PE2D.rotate_90(std_axes), (half_width_r2 - (half_width_r1 - d)/2) .* -i_cap)),
            (r2, r1, (half_width_r2 - d) .* -i_cap, std_axes, PE2D.Manifold(half_width_r1 + d, PE2D.rotate_90(std_axes), (half_width_r2 - (half_width_r1 + d)/2) .* -i_cap)),
            (r2, r1, (half_width_r2 - d) .* i_cap, std_axes, PE2D.Manifold(half_width_r1 + d, PE2D.rotate_minus_90(std_axes), (half_width_r2 - (half_width_r1 + d)/2) .* i_cap)),
            (r2, r1, (half_width_r2 + d) .* i_cap, std_axes, PE2D.Manifold(half_width_r1 - d, PE2D.rotate_minus_90(std_axes), (half_width_r2 - (half_width_r1 - d)/2) .* i_cap)),
            (r2, r1, (half_width_r2 + half_width_r1 - d) .* i_cap, std_axes, PE2D.Manifold(d, PE2D.rotate_minus_90(std_axes), (half_width_r2 - d/2) .* i_cap)),

            (r2, r1, top_right_r2 .- d, std_axes, PE2D.Manifold(half_height_r1 + d, std_axes, top_right_r2 .+ (half_width_r1 + d)/2 .* -i_cap .+ (half_height_r1 + d)/2 .* -j_cap)),
            (r2, r1, top_right_r2, std_axes, PE2D.Manifold(half_height_r1, std_axes, top_right_r2 .+ (half_width_r1/2) .* -i_cap .+ (half_height_r1/2) .* -j_cap)),
            (r2, r1, top_right_r2 .+ d, std_axes, PE2D.Manifold(half_height_r1 - d, std_axes, top_right_r2 .+ (half_width_r1 - d)/2 .* -i_cap .+ (half_height_r1 - d)/2 .* -j_cap)),

            # rotated_axes
            (r2, r1, (half_height_r2 - d) .* -j_cap .- PE2D.rotate(top_right_r1, rotated_axes), rotated_axes, PE2D.Manifold(d, PE2D.rotate_180(std_axes), ((zero(T) - d / tan(theta) + d * tan(theta)) ./ 3) .* i_cap .+ ((-half_height_r2 + d - half_height_r2 - half_height_r2) ./ 3) .* j_cap)),
            (r2, r1, (half_height_r2 - d) .* j_cap .+ PE2D.rotate(top_right_r1, rotated_axes), rotated_axes, PE2D.Manifold(d, std_axes, ((zero(T) + d / tan(theta) - d * tan(theta)) ./ 3) .* i_cap .+ ((half_height_r2 - d + half_height_r2 + half_height_r2) ./ 3) .* j_cap)),

            (r2, r1, (half_width_r2 + LA.norm(top_right_r1) * cos(-theta_r1 + theta) - d) .* -i_cap, rotated_axes, PE2D.Manifold(d, PE2D.rotate_90(std_axes), ((-half_width_r2 + d - half_width_r2 - half_width_r2) ./ 3) .* i_cap .+ (LA.norm(top_right_r1) * sin(-theta_r1 + theta) + (d / tan(theta) - d * tan(theta)) / 3) .* j_cap)),
            (r2, r1, (half_width_r2 + LA.norm(top_right_r1) * cos(-theta_r1 + theta) - d) .* i_cap, rotated_axes, PE2D.Manifold(d, PE2D.rotate_minus_90(std_axes), ((half_width_r2 - d + half_width_r2 + half_width_r2) ./ 3) .* i_cap .+ (LA.norm(top_right_r1) * sin(convert(T, pi - theta_r1) + theta) + (d * tan(theta) - d / tan(theta)) / 3) .* j_cap)),
            ]

            test_manifold_list(manifold_list)
        end
    end
end
