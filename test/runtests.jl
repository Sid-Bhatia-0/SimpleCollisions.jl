import LinearAlgebra as LA
import SimpleCollisions as SC
import StaticArrays as SA
import Test

function test_manifold_list_no_dir(manifold_list_no_dir)
    for (i, (a, b, pos_ba, value)) in enumerate(manifold_list_no_dir)
        manifold_ba = SC.Manifold(a, b, pos_ba)
        Test.@test SC.get_penetration(manifold_ba) ≈ SC.get_penetration(value)
        Test.@test SC.get_normal(manifold_ba) ≈ SC.get_normal(value)
        Test.@test SC.get_contact(manifold_ba) ≈ SC.get_contact(value)
    end
end

function test_manifold_list(manifold_list)
    for (i, (a, b, pos_ba, dir_ba, value)) in enumerate(manifold_list)
        manifold_ba = SC.Manifold(a, b, pos_ba, dir_ba)
        Test.@test SC.get_penetration(manifold_ba) ≈ SC.get_penetration(value)
        Test.@test SC.get_normal(manifold_ba) ≈ SC.get_normal(value)
        Test.@test SC.get_contact(manifold_ba) ≈ SC.get_contact(value)
    end
end

Test.@testset "SimpleCollisions.jl" begin
    T = Float32
    VecType = SC.Vector2D{T}

    origin = zero(VecType)
    d = convert(T, 0.01)

    std_dir = SC.Vector2D(one(T), zero(T))
    i_cap = SC.Vector2D(one(T), zero(T))
    j_cap = SC.Vector2D(zero(T), one(T))

    theta = convert(T, π / 6)
    rotated_dir = SC.Vector2D(cos(theta), sin(theta))
    theta_45 = convert(T, π / 4)
    unit_45 = (i_cap + j_cap) / convert(T, sqrt(2))

    point = SC.StdPoint{T}()

    half_length_l1 = one(T)
    l1 = SC.StdLine(half_length_l1)
    p1_l1 = SC.get_tail(l1)
    p2_l1 = SC.get_head(l1)

    half_length_l2 = convert(T, 2)
    l2 = SC.StdLine(half_length_l2)
    p1_l2 = SC.get_tail(l2)
    p2_l2 = SC.get_head(l2)

    r_c1 = one(T)
    c1 = SC.StdCircle(r_c1)

    r_c2 = convert(T, 2)
    c2 = SC.StdCircle(r_c2)

    half_width_r1 = one(T)
    half_height_r1 = convert(T, 0.5)
    r1 = SC.StdRect(half_width_r1, half_height_r1)
    top_right_r1 = SC.get_top_right(r1)
    theta_r1 = atan(half_height_r1, half_width_r1)

    half_width_r2 = convert(T, 2)
    half_height_r2 = one(T)
    r2 = SC.StdRect(half_width_r2, half_height_r2)
    top_right_r2 = SC.get_top_right(r2)
    theta_r2 = atan(half_height_r2, half_width_r2)

    Test.@testset "Area" begin
        Test.@testset "StdRect" begin
            Test.@test SC.get_area(r2) == convert(T, 8)
        end

        Test.@testset "StdCircle" begin
            Test.@test SC.get_area(c1) ≈ convert(T, π)
        end
    end

    Test.@testset "Collision detection" begin
        Test.@testset "StdLine vs. StdLine" begin
            Test.@testset "no dir" begin
                Test.@test SC.is_colliding(l1, l2, origin) == true
                Test.@test SC.is_colliding(l1, l2, (half_length_l1 + half_length_l2 + d) * -i_cap) == false
                Test.@test SC.is_colliding(l1, l2, (half_length_l1 + half_length_l2 - d) * -i_cap) == true
                Test.@test SC.is_colliding(l1, l2, (half_length_l1 + half_length_l2 - d) * i_cap) == true
                Test.@test SC.is_colliding(l1, l2, (half_length_l1 + half_length_l2 + d) * i_cap) == false
                Test.@test SC.is_colliding(l1, l2, d * -j_cap) == false
                Test.@test SC.is_colliding(l1, l2, d * j_cap) == false
            end

            Test.@testset "std_dir" begin
                Test.@test SC.is_colliding(l1, l2, origin, std_dir) == true
                Test.@test SC.is_colliding(l1, l2, (half_length_l1 + half_length_l2 + d) * -i_cap, std_dir) == false
                Test.@test SC.is_colliding(l1, l2, (half_length_l1 + half_length_l2 - d) * -i_cap, std_dir) == true
                Test.@test SC.is_colliding(l1, l2, (half_length_l1 + half_length_l2 - d) * i_cap, std_dir) == true
                Test.@test SC.is_colliding(l1, l2, (half_length_l1 + half_length_l2 + d) * i_cap, std_dir) == false
                Test.@test SC.is_colliding(l1, l2, d * -j_cap, std_dir) == false
                Test.@test SC.is_colliding(l1, l2, d * j_cap, std_dir) == false
            end

            Test.@testset "rotated_dir" begin
                Test.@test SC.is_colliding(l1, l2, origin, rotated_dir) == true
                Test.@test SC.is_colliding(l1, l2, (half_length_l1 + half_length_l2 + d) * -i_cap, rotated_dir) == false
                Test.@test SC.is_colliding(l1, l2, (half_length_l1 + half_length_l2 - d) * -i_cap, rotated_dir) == false
                Test.@test SC.is_colliding(l1, l2, (half_length_l1 + d) * -i_cap, rotated_dir) == false
                Test.@test SC.is_colliding(l1, l2, (half_length_l1 - d) * -i_cap, rotated_dir) == true
                Test.@test SC.is_colliding(l1, l2, (half_length_l1 - d) * i_cap, rotated_dir) == true
                Test.@test SC.is_colliding(l1, l2, (half_length_l1 + d) * i_cap, rotated_dir) == false
                Test.@test SC.is_colliding(l1, l2, (half_length_l1 + half_length_l2 - d) * i_cap, rotated_dir) == false
                Test.@test SC.is_colliding(l1, l2, (half_length_l1 + half_length_l2 + d) * i_cap, rotated_dir) == false
                Test.@test SC.is_colliding(l2, l1, (half_length_l1 * sin(theta) + d) * -j_cap, rotated_dir) == false
                Test.@test SC.is_colliding(l2, l1, (half_length_l1 * sin(theta) - d) * -j_cap, rotated_dir) == true
                Test.@test SC.is_colliding(l2, l1, d * -j_cap, rotated_dir) == true
                Test.@test SC.is_colliding(l2, l1, d * j_cap, rotated_dir) == true
                Test.@test SC.is_colliding(l2, l1, (half_length_l1 * sin(theta) - d) * j_cap, rotated_dir) == true
                Test.@test SC.is_colliding(l2, l1, (half_length_l1 * sin(theta) + d) * j_cap, rotated_dir) == false
            end
        end

        Test.@testset "StdCircle vs. StdPoint" begin
            Test.@testset "std_dir" begin
                Test.@test SC.is_colliding(c1, point, origin, std_dir) == true
                Test.@test SC.is_colliding(c1, point, (r_c1 + d) * -i_cap, std_dir) == false
                Test.@test SC.is_colliding(c1, point, (r_c1 - d) * -i_cap, std_dir) == true
                Test.@test SC.is_colliding(c1, point, (r_c1 - d) * i_cap, std_dir) == true
                Test.@test SC.is_colliding(c1, point, (r_c1 + d) * i_cap, std_dir) == false
                Test.@test SC.is_colliding(c1, point, (r_c1 + d) * -j_cap, std_dir) == false
                Test.@test SC.is_colliding(c1, point, (r_c1 - d) * -j_cap, std_dir) == true
                Test.@test SC.is_colliding(c1, point, (r_c1 - d) * j_cap, std_dir) == true
                Test.@test SC.is_colliding(c1, point, (r_c1 + d) * j_cap, std_dir) == false
                Test.@test SC.is_colliding(c1, point, (r_c1 + d) * unit_45, std_dir) == false
                Test.@test SC.is_colliding(c1, point, (r_c1 - d) * unit_45, std_dir) == true
            end

            Test.@testset "reverse check with std_dir" begin
                Test.@test SC.is_colliding(point, c1, origin, std_dir) == true
                Test.@test SC.is_colliding(point, c1, (r_c1 + d) * -i_cap, std_dir) == false
                Test.@test SC.is_colliding(point, c1, (r_c1 - d) * -i_cap, std_dir) == true
                Test.@test SC.is_colliding(point, c1, (r_c1 - d) * i_cap, std_dir) == true
                Test.@test SC.is_colliding(point, c1, (r_c1 + d) * i_cap, std_dir) == false
                Test.@test SC.is_colliding(point, c1, (r_c1 + d) * -j_cap, std_dir) == false
                Test.@test SC.is_colliding(point, c1, (r_c1 - d) * -j_cap, std_dir) == true
                Test.@test SC.is_colliding(point, c1, (r_c1 - d) * j_cap, std_dir) == true
                Test.@test SC.is_colliding(point, c1, (r_c1 + d) * j_cap, std_dir) == false
                Test.@test SC.is_colliding(point, c1, (r_c1 + d) * unit_45, std_dir) == false
                Test.@test SC.is_colliding(point, c1, (r_c1 - d) * unit_45, std_dir) == true
            end

            Test.@testset "no dir" begin
                Test.@test SC.is_colliding(c1, point, origin) == true
                Test.@test SC.is_colliding(c1, point, (r_c1 + d) * -i_cap) == false
                Test.@test SC.is_colliding(c1, point, (r_c1 - d) * -i_cap) == true
                Test.@test SC.is_colliding(c1, point, (r_c1 - d) * i_cap) == true
                Test.@test SC.is_colliding(c1, point, (r_c1 + d) * i_cap) == false
                Test.@test SC.is_colliding(c1, point, (r_c1 + d) * -j_cap) == false
                Test.@test SC.is_colliding(c1, point, (r_c1 - d) * -j_cap) == true
                Test.@test SC.is_colliding(c1, point, (r_c1 - d) * j_cap) == true
                Test.@test SC.is_colliding(c1, point, (r_c1 + d) * j_cap) == false
                Test.@test SC.is_colliding(c1, point, (r_c1 + d) * unit_45) == false
                Test.@test SC.is_colliding(c1, point, (r_c1 - d) * unit_45) == true
            end

            Test.@testset "reverse check with no dir" begin
                Test.@test SC.is_colliding(point, c1, origin) == true
                Test.@test SC.is_colliding(point, c1, (r_c1 + d) * -i_cap) == false
                Test.@test SC.is_colliding(point, c1, (r_c1 - d) * -i_cap) == true
                Test.@test SC.is_colliding(point, c1, (r_c1 - d) * i_cap) == true
                Test.@test SC.is_colliding(point, c1, (r_c1 + d) * i_cap) == false
                Test.@test SC.is_colliding(point, c1, (r_c1 + d) * -j_cap) == false
                Test.@test SC.is_colliding(point, c1, (r_c1 - d) * -j_cap) == true
                Test.@test SC.is_colliding(point, c1, (r_c1 - d) * j_cap) == true
                Test.@test SC.is_colliding(point, c1, (r_c1 + d) * j_cap) == false
                Test.@test SC.is_colliding(point, c1, (r_c1 + d) * unit_45) == false
                Test.@test SC.is_colliding(point, c1, (r_c1 - d) * unit_45) == true
            end
        end

        Test.@testset "StdCircle vs. StdLine" begin
            Test.@testset "std_dir" begin
                Test.@test SC.is_colliding(l1, c2, origin, std_dir) == true
                Test.@test SC.is_colliding(l2, c1, origin, std_dir) == true
                Test.@test SC.is_colliding(l2, c1, (half_length_l2 + r_c1 + d) * -i_cap, std_dir) == false
                Test.@test SC.is_colliding(l2, c1, (half_length_l2 + r_c1 - d) * -i_cap, std_dir) == true
                Test.@test SC.is_colliding(l2, c1, (half_length_l2 + r_c1 - d) * i_cap, std_dir) == true
                Test.@test SC.is_colliding(l2, c1, (half_length_l2 + r_c1 + d) * i_cap, std_dir) == false
                Test.@test SC.is_colliding(l2, c1, (r_c1 + d) * -j_cap, std_dir) == false
                Test.@test SC.is_colliding(l2, c1, (r_c1 - d) * -j_cap, std_dir) == true
                Test.@test SC.is_colliding(l2, c1, (r_c1 - d) * j_cap, std_dir) == true
                Test.@test SC.is_colliding(l2, c1, (r_c1 + d) * j_cap, std_dir) == false
            end

            Test.@testset "reverse check with std_dir" begin
                Test.@test SC.is_colliding(c2, l1, origin, std_dir) == true
                Test.@test SC.is_colliding(c1, l2, origin, std_dir) == true
                Test.@test SC.is_colliding(c1, l2, (half_length_l2 + r_c1 + d) * -i_cap, std_dir) == false
                Test.@test SC.is_colliding(c1, l2, (half_length_l2 + r_c1 - d) * -i_cap, std_dir) == true
                Test.@test SC.is_colliding(c1, l2, (half_length_l2 + r_c1 - d) * i_cap, std_dir) == true
                Test.@test SC.is_colliding(c1, l2, (half_length_l2 + r_c1 + d) * i_cap, std_dir) == false
                Test.@test SC.is_colliding(c1, l2, (r_c1 + d) * -j_cap, std_dir) == false
                Test.@test SC.is_colliding(c1, l2, (r_c1 - d) * -j_cap, std_dir) == true
                Test.@test SC.is_colliding(c1, l2, (r_c1 - d) * j_cap, std_dir) == true
                Test.@test SC.is_colliding(c1, l2, (r_c1 + d) * j_cap, std_dir) == false
            end

            Test.@testset "reverse check with rotated_dir" begin
                Test.@test SC.is_colliding(c2, l1, origin, rotated_dir) == true
                Test.@test SC.is_colliding(c1, l2, origin, rotated_dir) == true
                Test.@test SC.is_colliding(c1, l1, (sqrt(r_c1 ^ 2 - (half_length_l1 * sin(theta)) ^ 2) + half_length_l1 * cos(theta) + d) * -i_cap, rotated_dir) == false
                Test.@test SC.is_colliding(c1, l1, (sqrt(r_c1 ^ 2 - (half_length_l1 * sin(theta)) ^ 2) + half_length_l1 * cos(theta) - d) * -i_cap, rotated_dir) == true
                Test.@test SC.is_colliding(c1, l1, (sqrt(r_c1 ^ 2 - (half_length_l1 * sin(theta)) ^ 2) + half_length_l1 * cos(theta) - d) * i_cap, rotated_dir) == true
                Test.@test SC.is_colliding(c1, l1, (sqrt(r_c1 ^ 2 - (half_length_l1 * sin(theta)) ^ 2) + half_length_l1 * cos(theta) + d) * i_cap, rotated_dir) == false
                Test.@test SC.is_colliding(c1, l2, (r_c1 / cos(theta) + d) * -j_cap, rotated_dir) == false
                Test.@test SC.is_colliding(c1, l2, (r_c1 / cos(theta) - d) * -j_cap, rotated_dir) == true
                Test.@test SC.is_colliding(c1, l2, (r_c1 / cos(theta) - d) * j_cap, rotated_dir) == true
                Test.@test SC.is_colliding(c1, l2, (r_c1 / cos(theta) + d) * j_cap, rotated_dir) == false
            end

            Test.@testset "no dir" begin
                Test.@test SC.is_colliding(l1, c2, origin) == true
                Test.@test SC.is_colliding(l2, c1, origin) == true
                Test.@test SC.is_colliding(l2, c1, (half_length_l2 + r_c1 + d) * -i_cap) == false
                Test.@test SC.is_colliding(l2, c1, (half_length_l2 + r_c1 - d) * -i_cap) == true
                Test.@test SC.is_colliding(l2, c1, (half_length_l2 + r_c1 - d) * i_cap) == true
                Test.@test SC.is_colliding(l2, c1, (half_length_l2 + r_c1 + d) * i_cap) == false
                Test.@test SC.is_colliding(l2, c1, (r_c1 + d) * -j_cap) == false
                Test.@test SC.is_colliding(l2, c1, (r_c1 - d) * -j_cap) == true
                Test.@test SC.is_colliding(l2, c1, (r_c1 - d) * j_cap) == true
                Test.@test SC.is_colliding(l2, c1, (r_c1 + d) * j_cap) == false
            end

            Test.@testset "reverse check with no dir" begin
                Test.@test SC.is_colliding(c2, l1, origin) == true
                Test.@test SC.is_colliding(c1, l2, origin) == true
                Test.@test SC.is_colliding(c1, l2, (half_length_l2 + r_c1 + d) * -i_cap) == false
                Test.@test SC.is_colliding(c1, l2, (half_length_l2 + r_c1 - d) * -i_cap) == true
                Test.@test SC.is_colliding(c1, l2, (half_length_l2 + r_c1 - d) * i_cap) == true
                Test.@test SC.is_colliding(c1, l2, (half_length_l2 + r_c1 + d) * i_cap) == false
                Test.@test SC.is_colliding(c1, l2, (r_c1 + d) * -j_cap) == false
                Test.@test SC.is_colliding(c1, l2, (r_c1 - d) * -j_cap) == true
                Test.@test SC.is_colliding(c1, l2, (r_c1 - d) * j_cap) == true
                Test.@test SC.is_colliding(c1, l2, (r_c1 + d) * j_cap) == false
            end
        end

        Test.@testset "StdCircle vs. StdCircle" begin
            Test.@testset "std_dir" begin
                Test.@test SC.is_colliding(c1, c2, origin, std_dir) == true
                Test.@test SC.is_colliding(c1, c2, (r_c1 + r_c2 + d) * -i_cap, std_dir) == false
                Test.@test SC.is_colliding(c1, c2, (r_c1 + r_c2 - d) * -i_cap, std_dir) == true
                Test.@test SC.is_colliding(c1, c2, (r_c1 + r_c2 - d) * i_cap, std_dir) == true
                Test.@test SC.is_colliding(c1, c2, (r_c1 + r_c2 + d) * i_cap, std_dir) == false
                Test.@test SC.is_colliding(c1, c2, (r_c1 + r_c2 + d) * -j_cap, std_dir) == false
                Test.@test SC.is_colliding(c1, c2, (r_c1 + r_c2 - d) * -j_cap, std_dir) == true
                Test.@test SC.is_colliding(c1, c2, (r_c1 + r_c2 - d) * j_cap, std_dir) == true
                Test.@test SC.is_colliding(c1, c2, (r_c1 + r_c2 + d) * j_cap, std_dir) == false
                Test.@test SC.is_colliding(c1, c2, (r_c1 + r_c2 + d) * unit_45, std_dir) == false
                Test.@test SC.is_colliding(c1, c2, (r_c1 + r_c2 - d) * unit_45, std_dir) == true
            end

            Test.@testset "no dir" begin
                Test.@test SC.is_colliding(c1, c2, origin) == true
                Test.@test SC.is_colliding(c1, c2, (r_c1 + r_c2 + d) * -i_cap) == false
                Test.@test SC.is_colliding(c1, c2, (r_c1 + r_c2 - d) * -i_cap) == true
                Test.@test SC.is_colliding(c1, c2, (r_c1 + r_c2 - d) * i_cap) == true
                Test.@test SC.is_colliding(c1, c2, (r_c1 + r_c2 + d) * i_cap) == false
                Test.@test SC.is_colliding(c1, c2, (r_c1 + r_c2 + d) * -j_cap) == false
                Test.@test SC.is_colliding(c1, c2, (r_c1 + r_c2 - d) * -j_cap) == true
                Test.@test SC.is_colliding(c1, c2, (r_c1 + r_c2 - d) * j_cap) == true
                Test.@test SC.is_colliding(c1, c2, (r_c1 + r_c2 + d) * j_cap) == false
                Test.@test SC.is_colliding(c1, c2, (r_c1 + r_c2 + d) * unit_45) == false
                Test.@test SC.is_colliding(c1, c2, (r_c1 + r_c2 - d) * unit_45) == true
            end
        end

        Test.@testset "StdRect vs. StdPoint" begin
            Test.@testset "std_dir" begin
                Test.@test SC.is_colliding(r1, point, origin, std_dir) == true
                Test.@test SC.is_colliding(r1, point, (half_width_r1 + d) * -i_cap, std_dir) == false
                Test.@test SC.is_colliding(r1, point, (half_width_r1 - d) * -i_cap, std_dir) == true
                Test.@test SC.is_colliding(r1, point, (half_width_r1 - d) * i_cap, std_dir) == true
                Test.@test SC.is_colliding(r1, point, (half_width_r1 + d) * i_cap, std_dir) == false
                Test.@test SC.is_colliding(r1, point, (half_height_r1 + d) * -j_cap, std_dir) == false
                Test.@test SC.is_colliding(r1, point, (half_height_r1 - d) * -j_cap, std_dir) == true
                Test.@test SC.is_colliding(r1, point, (half_height_r1 - d) * j_cap, std_dir) == true
                Test.@test SC.is_colliding(r1, point, (half_height_r1 + d) * j_cap, std_dir) == false
                Test.@test SC.is_colliding(r1, point, top_right_r1 .+ d, std_dir) == false
                Test.@test SC.is_colliding(r1, point, top_right_r1 .- d, std_dir) == true
            end

            Test.@testset "reverse check with std_dir" begin
                Test.@test SC.is_colliding(point, r1, origin, std_dir) == true
                Test.@test SC.is_colliding(point, r1, (half_width_r1 + d) * -i_cap, std_dir) == false
                Test.@test SC.is_colliding(point, r1, (half_width_r1 - d) * -i_cap, std_dir) == true
                Test.@test SC.is_colliding(point, r1, (half_width_r1 - d) * i_cap, std_dir) == true
                Test.@test SC.is_colliding(point, r1, (half_width_r1 + d) * i_cap, std_dir) == false
                Test.@test SC.is_colliding(point, r1, (half_height_r1 + d) * -j_cap, std_dir) == false
                Test.@test SC.is_colliding(point, r1, (half_height_r1 - d) * -j_cap, std_dir) == true
                Test.@test SC.is_colliding(point, r1, (half_height_r1 - d) * j_cap, std_dir) == true
                Test.@test SC.is_colliding(point, r1, (half_height_r1 + d) * j_cap, std_dir) == false
                Test.@test SC.is_colliding(point, r1, top_right_r1 .+ d, std_dir) == false
                Test.@test SC.is_colliding(point, r1, top_right_r1 .- d, std_dir) == true
            end

            Test.@testset "reverse check with rotated_dir" begin
                Test.@test SC.is_colliding(point, r1, origin, rotated_dir) == true
                Test.@test SC.is_colliding(point, r1, (half_height_r1 / sin(theta) + d) * -i_cap, rotated_dir) == false
                Test.@test SC.is_colliding(point, r1, (half_height_r1 / sin(theta) - d) * -i_cap, rotated_dir) == true
                Test.@test SC.is_colliding(point, r1, (half_height_r1 / sin(theta) - d) * i_cap, rotated_dir) == true
                Test.@test SC.is_colliding(point, r1, (half_height_r1 / sin(theta) + d) * i_cap, rotated_dir) == false
                Test.@test SC.is_colliding(point, r1, (half_height_r1 / cos(theta) + d) * -j_cap, rotated_dir) == false
                Test.@test SC.is_colliding(point, r1, (half_height_r1 / cos(theta) - d) * -j_cap, rotated_dir) == true
                Test.@test SC.is_colliding(point, r1, (half_height_r1 / cos(theta) - d) * j_cap, rotated_dir) == true
                Test.@test SC.is_colliding(point, r1, (half_height_r1 / cos(theta) + d) * j_cap, rotated_dir) == false
            end

            Test.@testset "no dir" begin
                Test.@test SC.is_colliding(r1, point, origin) == true
                Test.@test SC.is_colliding(r1, point, (half_width_r1 + d) * -i_cap) == false
                Test.@test SC.is_colliding(r1, point, (half_width_r1 - d) * -i_cap) == true
                Test.@test SC.is_colliding(r1, point, (half_width_r1 - d) * i_cap) == true
                Test.@test SC.is_colliding(r1, point, (half_width_r1 + d) * i_cap) == false
                Test.@test SC.is_colliding(r1, point, (half_height_r1 + d) * -j_cap) == false
                Test.@test SC.is_colliding(r1, point, (half_height_r1 - d) * -j_cap) == true
                Test.@test SC.is_colliding(r1, point, (half_height_r1 - d) * j_cap) == true
                Test.@test SC.is_colliding(r1, point, (half_height_r1 + d) * j_cap) == false
                Test.@test SC.is_colliding(r1, point, top_right_r1 .+ d) == false
                Test.@test SC.is_colliding(r1, point, top_right_r1 .- d) == true
            end

            Test.@testset "reverse check with no dir" begin
                Test.@test SC.is_colliding(point, r1, origin) == true
                Test.@test SC.is_colliding(point, r1, (half_width_r1 + d) * -i_cap) == false
                Test.@test SC.is_colliding(point, r1, (half_width_r1 - d) * -i_cap) == true
                Test.@test SC.is_colliding(point, r1, (half_width_r1 - d) * i_cap) == true
                Test.@test SC.is_colliding(point, r1, (half_width_r1 + d) * i_cap) == false
                Test.@test SC.is_colliding(point, r1, (half_height_r1 + d) * -j_cap) == false
                Test.@test SC.is_colliding(point, r1, (half_height_r1 - d) * -j_cap) == true
                Test.@test SC.is_colliding(point, r1, (half_height_r1 - d) * j_cap) == true
                Test.@test SC.is_colliding(point, r1, (half_height_r1 + d) * j_cap) == false
                Test.@test SC.is_colliding(point, r1, top_right_r1 .+ d) == false
                Test.@test SC.is_colliding(point, r1, top_right_r1 .- d) == true
            end
        end

        Test.@testset "StdRect vs. StdLine" begin
            Test.@testset "std_dir" begin
                Test.@test SC.is_colliding(r1, l1, origin, std_dir) == true
                Test.@test SC.is_colliding(r1, l2, origin, std_dir) == true
                Test.@test SC.is_colliding(r1, l1, (half_width_r1 + half_length_l1 + d) * -i_cap, std_dir) == false
                Test.@test SC.is_colliding(r1, l1, (half_width_r1 + half_length_l1 - d) * -i_cap, std_dir) == true
                Test.@test SC.is_colliding(r1, l1, (half_width_r1 + half_length_l1 - d) * i_cap, std_dir) == true
                Test.@test SC.is_colliding(r1, l1, (half_width_r1 + half_length_l1 + d) * i_cap, std_dir) == false
                Test.@test SC.is_colliding(r1, l1, (half_height_r1 + d) * -j_cap, std_dir) == false
                Test.@test SC.is_colliding(r1, l1, (half_height_r1 - d) * -j_cap, std_dir) == true
                Test.@test SC.is_colliding(r1, l1, (half_height_r1 - d) * j_cap, std_dir) == true
                Test.@test SC.is_colliding(r1, l1, (half_height_r1 + d) * j_cap, std_dir) == false
            end

            Test.@testset "rotated_dir" begin
                Test.@test SC.is_colliding(r1, l1, origin, rotated_dir) == true
                Test.@test SC.is_colliding(r1, l2, origin, rotated_dir) == true
                Test.@test SC.is_colliding(r2, l1, (half_width_r2 + half_length_l1 * cos(theta) + d) * -i_cap, rotated_dir) == false
                Test.@test SC.is_colliding(r2, l1, (half_width_r2 + half_length_l1 * cos(theta) - d) * -i_cap, rotated_dir) == true
                Test.@test SC.is_colliding(r2, l1, (half_width_r2 + half_length_l1 * cos(theta) - d) * i_cap, rotated_dir) == true
                Test.@test SC.is_colliding(r2, l1, (half_width_r2 + half_length_l1 * cos(theta) + d) * i_cap, rotated_dir) == false
                Test.@test SC.is_colliding(r2, l1, (half_height_r2 + half_length_l1 * sin(theta) + d) * -j_cap, rotated_dir) == false
                Test.@test SC.is_colliding(r2, l1, (half_height_r2 + half_length_l1 * sin(theta) - d) * -j_cap, rotated_dir) == true
                Test.@test SC.is_colliding(r2, l1, (half_height_r2 + half_length_l1 * sin(theta) - d) * j_cap, rotated_dir) == true
                Test.@test SC.is_colliding(r2, l1, (half_height_r2 + half_length_l1 * sin(theta) + d) * j_cap, rotated_dir) == false
            end

            Test.@testset "reverse check with std_dir" begin
                Test.@test SC.is_colliding(l1, r1, origin, std_dir) == true
                Test.@test SC.is_colliding(l2, r1, origin, std_dir) == true
                Test.@test SC.is_colliding(l1, r1, (half_width_r1 + half_length_l1 + d) * -i_cap, std_dir) == false
                Test.@test SC.is_colliding(l1, r1, (half_width_r1 + half_length_l1 - d) * -i_cap, std_dir) == true
                Test.@test SC.is_colliding(l1, r1, (half_width_r1 + half_length_l1 - d) * i_cap, std_dir) == true
                Test.@test SC.is_colliding(l1, r1, (half_width_r1 + half_length_l1 + d) * i_cap, std_dir) == false
                Test.@test SC.is_colliding(l1, r1, (half_height_r1 + d) * -j_cap, std_dir) == false
                Test.@test SC.is_colliding(l1, r1, (half_height_r1 - d) * -j_cap, std_dir) == true
                Test.@test SC.is_colliding(l1, r1, (half_height_r1 - d) * j_cap, std_dir) == true
                Test.@test SC.is_colliding(l1, r1, (half_height_r1 + d) * j_cap, std_dir) == false
            end

            Test.@testset "reverse check with rotated_dir" begin
                Test.@test SC.is_colliding(l1, r1, origin, rotated_dir) == true
                Test.@test SC.is_colliding(l2, r1, origin, rotated_dir) == true
                Test.@test SC.is_colliding(l2, r2, (half_height_r2 / sin(theta) + half_length_l2 + d) * -i_cap, rotated_dir) == false
                Test.@test SC.is_colliding(l2, r2, (half_height_r2 / sin(theta) + half_length_l2 - d) * -i_cap, rotated_dir) == true
                Test.@test SC.is_colliding(l2, r2, (half_height_r2 / sin(theta) + half_length_l2 - d) * i_cap, rotated_dir) == true
                Test.@test SC.is_colliding(l2, r2, (half_height_r2 / sin(theta) + half_length_l2 + d) * i_cap, rotated_dir) == false
                Test.@test SC.is_colliding(l2, r2, (half_width_r2 * sin(theta) + half_height_r2 * cos(theta) + d) * -j_cap, rotated_dir) == false
                Test.@test SC.is_colliding(l2, r2, (half_width_r2 * sin(theta) + half_height_r2 * cos(theta) - d) * -j_cap, rotated_dir) == true
                Test.@test SC.is_colliding(l2, r2, (half_width_r2 * sin(theta) + half_height_r2 * cos(theta) - d) * j_cap, rotated_dir) == true
                Test.@test SC.is_colliding(l2, r2, (half_width_r2 * sin(theta) + half_height_r2 * cos(theta) + d) * j_cap, rotated_dir) == false
            end

            Test.@testset "no dir" begin
                Test.@test SC.is_colliding(r1, l1, origin) == true
                Test.@test SC.is_colliding(r1, l2, origin) == true
                Test.@test SC.is_colliding(r1, l1, (half_width_r1 + half_length_l1 + d) * -i_cap) == false
                Test.@test SC.is_colliding(r1, l1, (half_width_r1 + half_length_l1 - d) * -i_cap) == true
                Test.@test SC.is_colliding(r1, l1, (half_width_r1 + half_length_l1 - d) * i_cap) == true
                Test.@test SC.is_colliding(r1, l1, (half_width_r1 + half_length_l1 + d) * i_cap) == false
                Test.@test SC.is_colliding(r1, l1, (half_height_r1 + d) * -j_cap) == false
                Test.@test SC.is_colliding(r1, l1, (half_height_r1 - d) * -j_cap) == true
                Test.@test SC.is_colliding(r1, l1, (half_height_r1 - d) * j_cap) == true
                Test.@test SC.is_colliding(r1, l1, (half_height_r1 + d) * j_cap) == false
            end

            Test.@testset "reverse check with no dir" begin
                Test.@test SC.is_colliding(l1, r1, origin) == true
                Test.@test SC.is_colliding(l2, r1, origin) == true
                Test.@test SC.is_colliding(l1, r1, (half_width_r1 + half_length_l1 + d) * -i_cap) == false
                Test.@test SC.is_colliding(l1, r1, (half_width_r1 + half_length_l1 - d) * -i_cap) == true
                Test.@test SC.is_colliding(l1, r1, (half_width_r1 + half_length_l1 - d) * i_cap) == true
                Test.@test SC.is_colliding(l1, r1, (half_width_r1 + half_length_l1 + d) * i_cap) == false
                Test.@test SC.is_colliding(l1, r1, (half_height_r1 + d) * -j_cap) == false
                Test.@test SC.is_colliding(l1, r1, (half_height_r1 - d) * -j_cap) == true
                Test.@test SC.is_colliding(l1, r1, (half_height_r1 - d) * j_cap) == true
                Test.@test SC.is_colliding(l1, r1, (half_height_r1 + d) * j_cap) == false
            end
        end

        Test.@testset "StdRect vs. StdCircle" begin
            Test.@testset "std_dir" begin
                Test.@test SC.is_colliding(r1, c1, origin, std_dir) == true
                Test.@test SC.is_colliding(r1, c1, (half_width_r1 + r_c1 + d) * -i_cap, std_dir) == false
                Test.@test SC.is_colliding(r1, c1, (half_width_r1 + r_c1 - d) * -i_cap, std_dir) == true
                Test.@test SC.is_colliding(r1, c1, (half_width_r1 + r_c1 - d) * i_cap, std_dir) == true
                Test.@test SC.is_colliding(r1, c1, (half_width_r1 + r_c1 + d) * i_cap, std_dir) == false
                Test.@test SC.is_colliding(r1, c1, (half_height_r1 + r_c1 + d) * -j_cap, std_dir) == false
                Test.@test SC.is_colliding(r1, c1, (half_height_r1 + r_c1 - d) * -j_cap, std_dir) == true
                Test.@test SC.is_colliding(r1, c1, (half_height_r1 + r_c1 - d) * j_cap, std_dir) == true
                Test.@test SC.is_colliding(r1, c1, (half_height_r1 + r_c1 + d) * -j_cap, std_dir) == false
                Test.@test SC.is_colliding(r1, c1, top_right_r1 + (r_c1 + d) * unit_45, std_dir) == false
                Test.@test SC.is_colliding(r1, c1, top_right_r1 + (r_c1 - d) * unit_45, std_dir) == true
            end

            Test.@testset "reverse check with std_dir" begin
                Test.@test SC.is_colliding(c1, r1, origin, std_dir) == true
                Test.@test SC.is_colliding(c1, r1, (half_width_r1 + r_c1 + d) * -i_cap, std_dir) == false
                Test.@test SC.is_colliding(c1, r1, (half_width_r1 + r_c1 - d) * -i_cap, std_dir) == true
                Test.@test SC.is_colliding(c1, r1, (half_width_r1 + r_c1 - d) * i_cap, std_dir) == true
                Test.@test SC.is_colliding(c1, r1, (half_width_r1 + r_c1 + d) * i_cap, std_dir) == false
                Test.@test SC.is_colliding(c1, r1, (half_height_r1 + r_c1 + d) * -j_cap, std_dir) == false
                Test.@test SC.is_colliding(c1, r1, (half_height_r1 + r_c1 - d) * -j_cap, std_dir) == true
                Test.@test SC.is_colliding(c1, r1, (half_height_r1 + r_c1 - d) * j_cap, std_dir) == true
                Test.@test SC.is_colliding(c1, r1, (half_height_r1 + r_c1 + d) * -j_cap, std_dir) == false
                Test.@test SC.is_colliding(c1, r1, top_right_r1 + (r_c1 + d) * unit_45, std_dir) == false
                Test.@test SC.is_colliding(c1, r1, top_right_r1 + (r_c1 - d) * unit_45, std_dir) == true
            end

            Test.@testset "reverse check with rotated_dir" begin
                Test.@test SC.is_colliding(c2, r2, origin, rotated_dir) == true
                Test.@test SC.is_colliding(c2, r2, (sqrt(r_c2 ^ 2 - (half_width_r2 * sin(theta) - half_height_r2 * cos(theta)) ^ 2) + half_width_r2 * cos(theta) + half_height_r2 * sin(theta) + d) * -i_cap, rotated_dir) == false
                Test.@test SC.is_colliding(c2, r2, (sqrt(r_c2 ^ 2 - (half_width_r2 * sin(theta) - half_height_r2 * cos(theta)) ^ 2) + half_width_r2 * cos(theta) + half_height_r2 * sin(theta) - d) * -i_cap, rotated_dir) == true
                Test.@test SC.is_colliding(c2, r2, (sqrt(r_c2 ^ 2 - (half_width_r2 * sin(theta) - half_height_r2 * cos(theta)) ^ 2) + half_width_r2 * cos(theta) + half_height_r2 * sin(theta) - d) * i_cap, rotated_dir) == true
                Test.@test SC.is_colliding(c2, r2, (sqrt(r_c2 ^ 2 - (half_width_r2 * sin(theta) - half_height_r2 * cos(theta)) ^ 2) + half_width_r2 * cos(theta) + half_height_r2 * sin(theta) + d) * i_cap, rotated_dir) == false
                Test.@test SC.is_colliding(c2, r2, ((r_c2 + half_height_r2) / cos(theta) + d) * -j_cap, rotated_dir) == false
                Test.@test SC.is_colliding(c2, r2, ((r_c2 + half_height_r2) / cos(theta) - d) * -j_cap, rotated_dir) == true
                Test.@test SC.is_colliding(c2, r2, ((r_c2 + half_height_r2) / cos(theta) - d) * j_cap, rotated_dir) == true
                Test.@test SC.is_colliding(c2, r2, ((r_c2 + half_height_r2) / cos(theta) + d) * j_cap, rotated_dir) == false
            end

            Test.@testset "no dir" begin
                Test.@test SC.is_colliding(r1, c1, origin) == true
                Test.@test SC.is_colliding(r1, c1, (half_width_r1 + r_c1 + d) * -i_cap) == false
                Test.@test SC.is_colliding(r1, c1, (half_width_r1 + r_c1 - d) * -i_cap) == true
                Test.@test SC.is_colliding(r1, c1, (half_width_r1 + r_c1 - d) * i_cap) == true
                Test.@test SC.is_colliding(r1, c1, (half_width_r1 + r_c1 + d) * i_cap) == false
                Test.@test SC.is_colliding(r1, c1, (half_height_r1 + r_c1 + d) * -j_cap) == false
                Test.@test SC.is_colliding(r1, c1, (half_height_r1 + r_c1 - d) * -j_cap) == true
                Test.@test SC.is_colliding(r1, c1, (half_height_r1 + r_c1 - d) * j_cap) == true
                Test.@test SC.is_colliding(r1, c1, (half_height_r1 + r_c1 + d) * -j_cap) == false
                Test.@test SC.is_colliding(r1, c1, top_right_r1 + (r_c1 + d) * unit_45) == false
                Test.@test SC.is_colliding(r1, c1, top_right_r1 + (r_c1 - d) * unit_45) == true
            end

            Test.@testset "reverse check with no dir" begin
                Test.@test SC.is_colliding(c1, r1, origin) == true
                Test.@test SC.is_colliding(c1, r1, (half_width_r1 + r_c1 + d) * -i_cap) == false
                Test.@test SC.is_colliding(c1, r1, (half_width_r1 + r_c1 - d) * -i_cap) == true
                Test.@test SC.is_colliding(c1, r1, (half_width_r1 + r_c1 - d) * i_cap) == true
                Test.@test SC.is_colliding(c1, r1, (half_width_r1 + r_c1 + d) * i_cap) == false
                Test.@test SC.is_colliding(c1, r1, (half_height_r1 + r_c1 + d) * -j_cap) == false
                Test.@test SC.is_colliding(c1, r1, (half_height_r1 + r_c1 - d) * -j_cap) == true
                Test.@test SC.is_colliding(c1, r1, (half_height_r1 + r_c1 - d) * j_cap) == true
                Test.@test SC.is_colliding(c1, r1, (half_height_r1 + r_c1 + d) * -j_cap) == false
                Test.@test SC.is_colliding(c1, r1, top_right_r1 + (r_c1 + d) * unit_45) == false
                Test.@test SC.is_colliding(c1, r1, top_right_r1 + (r_c1 - d) * unit_45) == true
            end
        end

        Test.@testset "Rect2D vs. Rect2D" begin
            Test.@testset "std_dir" begin
                Test.@test SC.is_colliding(r2, r1, origin, std_dir) == true
                Test.@test SC.is_colliding(r2, r1, (half_width_r1 + half_width_r2 + d) * -i_cap, std_dir) == false
                Test.@test SC.is_colliding(r2, r1, (half_width_r1 + half_width_r2 - d) * -i_cap, std_dir) == true
                Test.@test SC.is_colliding(r2, r1, (half_width_r1 + half_width_r2 - d) * i_cap, std_dir) == true
                Test.@test SC.is_colliding(r2, r1, (half_width_r1 + half_width_r2 + d) * i_cap, std_dir) == false
                Test.@test SC.is_colliding(r2, r1, (half_height_r1 + half_height_r2 + d) * -j_cap, std_dir) == false
                Test.@test SC.is_colliding(r2, r1, (half_height_r1 + half_height_r2 - d) * -j_cap, std_dir) == true
                Test.@test SC.is_colliding(r2, r1, (half_height_r1 + half_height_r2 - d) * j_cap, std_dir) == true
                Test.@test SC.is_colliding(r2, r1, (half_height_r1 + half_height_r2 + d) * j_cap, std_dir) == false
                Test.@test SC.is_colliding(r2, r1, top_right_r1 + top_right_r2 .+ d, std_dir) == false
                Test.@test SC.is_colliding(r2, r1, top_right_r1 + top_right_r2 .- d, std_dir) == true
            end

            Test.@testset "rotated_dir" begin
                Test.@test SC.is_colliding(r2, r1, origin, rotated_dir) == true
                Test.@test SC.is_colliding(r2, r1, (half_width_r2 + half_width_r1 * cos(theta) + half_height_r1 * sin(theta) + d) * -i_cap, rotated_dir) == false
                Test.@test SC.is_colliding(r2, r1, (half_width_r2 + half_width_r1 * cos(theta) + half_height_r1 * sin(theta) - d) * -i_cap, rotated_dir) == true
                Test.@test SC.is_colliding(r2, r1, (half_width_r2 + half_width_r1 * cos(theta) + half_height_r1 * sin(theta) - d) * i_cap, rotated_dir) == true
                Test.@test SC.is_colliding(r2, r1, (half_width_r2 + half_width_r1 * cos(theta) + half_height_r1 * sin(theta) + d) * i_cap, rotated_dir) == false
                Test.@test SC.is_colliding(r2, r1, (half_height_r2 + half_width_r1 * sin(theta) + half_height_r1 * cos(theta) + d) * -j_cap, rotated_dir) == false
                Test.@test SC.is_colliding(r2, r1, (half_height_r2 + half_width_r1 * sin(theta) + half_height_r1 * cos(theta) - d) * -j_cap, rotated_dir) == true
                Test.@test SC.is_colliding(r2, r1, (half_height_r2 + half_width_r1 * sin(theta) + half_height_r1 * cos(theta) - d) * j_cap, rotated_dir) == true
                Test.@test SC.is_colliding(r2, r1, (half_height_r2 + half_width_r1 * sin(theta) + half_height_r1 * cos(theta) + d) * j_cap, rotated_dir) == false
            end

            # collision_list = [
            # # std_dir
            # (r2, r1, origin, std_dir, true),

            # (r2, r1, (half_width_r1 + half_width_r2 + d) * -i_cap, std_dir, false),
            # (r2, r1, (half_width_r1 + half_width_r2 - d) * -i_cap, std_dir, true),
            # (r2, r1, (half_width_r1 + half_width_r2 - d) * i_cap, std_dir, true),
            # (r2, r1, (half_width_r1 + half_width_r2 + d) * i_cap, std_dir, false),

            # (r2, r1, (half_height_r1 + half_height_r2 + d) * -j_cap, std_dir, false),
            # (r2, r1, (half_height_r1 + half_height_r2 - d) * -j_cap, std_dir, true),
            # (r2, r1, (half_height_r1 + half_height_r2 - d) * j_cap, std_dir, true),
            # (r2, r1, (half_height_r1 + half_height_r2 + d) * j_cap, std_dir, false),

            # (r2, r1, top_right_r1 + top_right_r2 .+ d, std_dir, false),
            # (r2, r1, top_right_r1 + top_right_r2 .- d, std_dir, true),

            # # rotated_dir
            # (r2, r1, origin, rotated_dir, true),

            # (r2, r1, (half_width_r2 + half_width_r1 * cos(theta) + half_height_r1 * sin(theta) + d) * -i_cap, rotated_dir, false),
            # (r2, r1, (half_width_r2 + half_width_r1 * cos(theta) + half_height_r1 * sin(theta) - d) * -i_cap, rotated_dir, true),
            # (r2, r1, (half_width_r2 + half_width_r1 * cos(theta) + half_height_r1 * sin(theta) - d) * i_cap, rotated_dir, true),
            # (r2, r1, (half_width_r2 + half_width_r1 * cos(theta) + half_height_r1 * sin(theta) + d) * i_cap, rotated_dir, false),

            # (r2, r1, (half_height_r2 + half_width_r1 * sin(theta) + half_height_r1 * cos(theta) + d) * -j_cap, rotated_dir, false),
            # (r2, r1, (half_height_r2 + half_width_r1 * sin(theta) + half_height_r1 * cos(theta) - d) * -j_cap, rotated_dir, true),
            # (r2, r1, (half_height_r2 + half_width_r1 * sin(theta) + half_height_r1 * cos(theta) - d) * j_cap, rotated_dir, true),
            # (r2, r1, (half_height_r2 + half_width_r1 * sin(theta) + half_height_r1 * cos(theta) + d) * j_cap, rotated_dir, false),
            # ]

            # test_collision_list(collision_list)

            Test.@testset "no dir" begin
                Test.@test SC.is_colliding(r2, r1, origin) == true
                Test.@test SC.is_colliding(r2, r1, (half_width_r1 + half_width_r2 + d) * -i_cap) == false
                Test.@test SC.is_colliding(r2, r1, (half_width_r1 + half_width_r2 - d) * -i_cap) == true
                Test.@test SC.is_colliding(r2, r1, (half_width_r1 + half_width_r2 - d) * i_cap) == true
                Test.@test SC.is_colliding(r2, r1, (half_width_r1 + half_width_r2 + d) * i_cap) == false
                Test.@test SC.is_colliding(r2, r1, (half_height_r1 + half_height_r2 + d) * -j_cap) == false
                Test.@test SC.is_colliding(r2, r1, (half_height_r1 + half_height_r2 - d) * -j_cap) == true
                Test.@test SC.is_colliding(r2, r1, (half_height_r1 + half_height_r2 - d) * j_cap) == true
                Test.@test SC.is_colliding(r2, r1, (half_height_r1 + half_height_r2 + d) * j_cap) == false
                Test.@test SC.is_colliding(r2, r1, top_right_r1 + top_right_r2 .+ d) == false
                Test.@test SC.is_colliding(r2, r1, top_right_r1 + top_right_r2 .- d) == true
            end

            # collision_list_no_dir = [
            # # std_dir
            # (r2, r1, origin, true),

            # (r2, r1, (half_width_r1 + half_width_r2 + d) * -i_cap, false),
            # (r2, r1, (half_width_r1 + half_width_r2 - d) * -i_cap, true),
            # (r2, r1, (half_width_r1 + half_width_r2 - d) * i_cap, true),
            # (r2, r1, (half_width_r1 + half_width_r2 + d) * i_cap, false),

            # (r2, r1, (half_height_r1 + half_height_r2 + d) * -j_cap, false),
            # (r2, r1, (half_height_r1 + half_height_r2 - d) * -j_cap, true),
            # (r2, r1, (half_height_r1 + half_height_r2 - d) * j_cap, true),
            # (r2, r1, (half_height_r1 + half_height_r2 + d) * j_cap, false),

            # (r2, r1, top_right_r1 + top_right_r2 .+ d, false),
            # (r2, r1, top_right_r1 + top_right_r2 .- d, true),
            # ]

            # test_collision_list_no_dir(collision_list_no_dir)
        end
    end

    Test.@testset "Manifold generation" begin
        Test.@testset "StdCircle vs. StdCircle" begin
            manifold_list = [
            # std_dir
            (c1, c2, (r_c1 + r_c2 - d) * -i_cap, std_dir, SC.Manifold(d, -i_cap, (r_c1 - d / 2) * -i_cap)),
            (c1, c2, r_c2 * -i_cap, std_dir, SC.Manifold(r_c1, -i_cap, (r_c1 / 2) * -i_cap)),
            (c1, c2, r_c2 * i_cap, std_dir, SC.Manifold(r_c1, i_cap, (r_c1 / 2) * i_cap)),
            (c1, c2, (r_c1 + r_c2 - d) * i_cap, std_dir, SC.Manifold(d, i_cap, (r_c1 - d / 2) * i_cap)),

            (c1, c2, (r_c1 + r_c2 - d) * -j_cap, std_dir, SC.Manifold(d, -j_cap, (r_c1 - d / 2) * -j_cap)),
            (c1, c2, r_c2 * -j_cap, std_dir, SC.Manifold(r_c1, -j_cap, (r_c1 / 2) * -j_cap)),
            (c1, c2, r_c2 * j_cap, std_dir, SC.Manifold(r_c1, j_cap, (r_c1 / 2) * j_cap)),
            (c1, c2, (r_c1 + r_c2 - d) * j_cap, std_dir, SC.Manifold(d, j_cap, (r_c1 - d / 2) * j_cap)),

            (c1, c2, (r_c1 + r_c2 - d) * unit_45, std_dir, SC.Manifold(d, unit_45, (r_c1 - d / 2) * unit_45)),
            (c1, c2, r_c2 * unit_45, std_dir, SC.Manifold(r_c1, unit_45, (r_c1 / 2) * unit_45)),
            ]

            test_manifold_list(manifold_list)
        end

        Test.@testset "StdRect vs. StdCircle" begin
            manifold_list_no_dir = [
            # std_dir
            (r1, c1, (half_width_r1 + r_c1 - d) * -i_cap, SC.Manifold(d, -i_cap, (half_width_r1 - d / 2) * -i_cap)),
            (r1, c1, (half_width_r1 + d) * -i_cap, SC.Manifold(r_c1 - d, -i_cap, (half_width_r1 - (r_c1 - d) / 2) * -i_cap)),
            (r1, c1, (half_width_r1 - d) * -i_cap, SC.Manifold(r_c1 + d, -i_cap, (half_width_r1 - (r_c1 + d) / 2) * -i_cap)),
            (r1, c1, (half_width_r1 - d) * i_cap, SC.Manifold(r_c1 + d, i_cap, (half_width_r1 - (r_c1 + d) / 2) * i_cap)),
            (r1, c1, (half_width_r1 + d) * i_cap, SC.Manifold(r_c1 - d, i_cap, (half_width_r1 - (r_c1 - d) / 2) * i_cap)),
            (r1, c1, (half_width_r1 + r_c1 - d) * i_cap, SC.Manifold(d, i_cap, (half_width_r1 - d / 2) * i_cap)),

            (r1, c1, (half_height_r1 + r_c1 - d) * -j_cap, SC.Manifold(d, -j_cap, (half_height_r1 - d / 2) * -j_cap)),
            (r1, c1, (half_height_r1 + d) * -j_cap, SC.Manifold(r_c1 - d, -j_cap, (half_height_r1 - (r_c1 - d) / 2) * -j_cap)),
            (r1, c1, (half_height_r1 - d) * -j_cap, SC.Manifold(r_c1 + d, -j_cap, (half_height_r1 - (r_c1 + d) / 2) * -j_cap)),
            (r1, c1, (half_height_r1 - d) * j_cap, SC.Manifold(r_c1 + d, j_cap, (half_height_r1 - (r_c1 + d) / 2) * j_cap)),
            (r1, c1, (half_height_r1 + d) * j_cap, SC.Manifold(r_c1 - d, j_cap, (half_height_r1 - (r_c1 - d) / 2) * j_cap)),
            (r1, c1, (half_height_r1 + r_c1 - d) * j_cap, SC.Manifold(d, j_cap, (half_height_r1 - d / 2) * j_cap)),

            (r1, c1, top_right_r1 + (r_c1 - d) * unit_45, SC.Manifold(d, unit_45, top_right_r1 + (d / 2) * -unit_45)),

            # reverse check with std_dir
            (c1, r1, (half_width_r1 + r_c1 - d) * -i_cap, SC.Manifold(d, -i_cap, (r_c1 - d / 2) * -i_cap)),
            (c1, r1, (half_width_r1 + d) * -i_cap, SC.Manifold(r_c1 - d, -i_cap, (r_c1 - (r_c1 - d) / 2) * -i_cap)),
            (c1, r1, (half_width_r1 - d) * -i_cap, SC.Manifold(r_c1 + d, -i_cap, (r_c1 - (r_c1 + d) / 2) * -i_cap)),
            (c1, r1, (half_width_r1 - d) * i_cap, SC.Manifold(r_c1 + d, i_cap, (r_c1 - (r_c1 + d) / 2) * i_cap)),
            (c1, r1, (half_width_r1 + d) * i_cap, SC.Manifold(r_c1 - d, i_cap, (r_c1 - (r_c1 - d) / 2) * i_cap)),
            (c1, r1, (half_width_r1 + r_c1 - d) * i_cap, SC.Manifold(d, i_cap, (r_c1 - d / 2) * i_cap)),

            (c1, r1, (half_height_r1 + r_c1 - d) * -j_cap, SC.Manifold(d, -j_cap, (r_c1 - d / 2) * -j_cap)),
            (c1, r1, (half_height_r1 + d) * -j_cap, SC.Manifold(r_c1 - d, -j_cap, (r_c1 - (r_c1 - d) / 2) * -j_cap)),
            (c1, r1, (half_height_r1 - d) * -j_cap, SC.Manifold(r_c1 + d, -j_cap, (r_c1 - (r_c1 + d) / 2) * -j_cap)),
            (c1, r1, (half_height_r1 - d) * j_cap, SC.Manifold(r_c1 + d, j_cap, (r_c1 - (r_c1 + d) / 2) * j_cap)),
            (c1, r1, (half_height_r1 + d) * j_cap, SC.Manifold(r_c1 - d, j_cap, (r_c1 - (r_c1 - d) / 2) * j_cap)),
            (c1, r1, (half_height_r1 + r_c1 - d) * j_cap, SC.Manifold(d, j_cap, (r_c1 - d / 2) * j_cap)),

            (c1, r1, top_right_r1 + (r_c1 - d) * unit_45, SC.Manifold(d, unit_45, (r_c1 - d / 2) * unit_45)),
            ]

            test_manifold_list_no_dir(manifold_list_no_dir)

            manifold_list = [
            # std_dir
            (r1, c1, (half_width_r1 + r_c1 - d) * -i_cap, std_dir, SC.Manifold(d, -i_cap, (half_width_r1 - d / 2) * -i_cap)),
            (r1, c1, (half_width_r1 + d) * -i_cap, std_dir, SC.Manifold(r_c1 - d, -i_cap, (half_width_r1 - (r_c1 - d) / 2) * -i_cap)),
            (r1, c1, (half_width_r1 - d) * -i_cap, std_dir, SC.Manifold(r_c1 + d, -i_cap, (half_width_r1 - (r_c1 + d) / 2) * -i_cap)),
            (r1, c1, (half_width_r1 - d) * i_cap, std_dir, SC.Manifold(r_c1 + d, i_cap, (half_width_r1 - (r_c1 + d) / 2) * i_cap)),
            (r1, c1, (half_width_r1 + d) * i_cap, std_dir, SC.Manifold(r_c1 - d, i_cap, (half_width_r1 - (r_c1 - d) / 2) * i_cap)),
            (r1, c1, (half_width_r1 + r_c1 - d) * i_cap, std_dir, SC.Manifold(d, i_cap, (half_width_r1 - d / 2) * i_cap)),

            (r1, c1, (half_height_r1 + r_c1 - d) * -j_cap, std_dir, SC.Manifold(d, -j_cap, (half_height_r1 - d / 2) * -j_cap)),
            (r1, c1, (half_height_r1 + d) * -j_cap, std_dir, SC.Manifold(r_c1 - d, -j_cap, (half_height_r1 - (r_c1 - d) / 2) * -j_cap)),
            (r1, c1, (half_height_r1 - d) * -j_cap, std_dir, SC.Manifold(r_c1 + d, -j_cap, (half_height_r1 - (r_c1 + d) / 2) * -j_cap)),
            (r1, c1, (half_height_r1 - d) * j_cap, std_dir, SC.Manifold(r_c1 + d, j_cap, (half_height_r1 - (r_c1 + d) / 2) * j_cap)),
            (r1, c1, (half_height_r1 + d) * j_cap, std_dir, SC.Manifold(r_c1 - d, j_cap, (half_height_r1 - (r_c1 - d) / 2) * j_cap)),
            (r1, c1, (half_height_r1 + r_c1 - d) * j_cap, std_dir, SC.Manifold(d, j_cap, (half_height_r1 - d / 2) * j_cap)),

            (r1, c1, top_right_r1 + (r_c1 - d) * unit_45, std_dir, SC.Manifold(d, unit_45, top_right_r1 + (d / 2) * -unit_45)),

            # reverse check with std_dir
            (c1, r1, (half_width_r1 + r_c1 - d) * -i_cap, std_dir, SC.Manifold(d, -i_cap, (r_c1 - d / 2) * -i_cap)),
            (c1, r1, (half_width_r1 + d) * -i_cap, std_dir, SC.Manifold(r_c1 - d, -i_cap, (r_c1 - (r_c1 - d) / 2) * -i_cap)),
            (c1, r1, (half_width_r1 - d) * -i_cap, std_dir, SC.Manifold(r_c1 + d, -i_cap, (r_c1 - (r_c1 + d) / 2) * -i_cap)),
            (c1, r1, (half_width_r1 - d) * i_cap, std_dir, SC.Manifold(r_c1 + d, i_cap, (r_c1 - (r_c1 + d) / 2) * i_cap)),
            (c1, r1, (half_width_r1 + d) * i_cap, std_dir, SC.Manifold(r_c1 - d, i_cap, (r_c1 - (r_c1 - d) / 2) * i_cap)),
            (c1, r1, (half_width_r1 + r_c1 - d) * i_cap, std_dir, SC.Manifold(d, i_cap, (r_c1 - d / 2) * i_cap)),

            (c1, r1, (half_height_r1 + r_c1 - d) * -j_cap, std_dir, SC.Manifold(d, -j_cap, (r_c1 - d / 2) * -j_cap)),
            (c1, r1, (half_height_r1 + d) * -j_cap, std_dir, SC.Manifold(r_c1 - d, -j_cap, (r_c1 - (r_c1 - d) / 2) * -j_cap)),
            (c1, r1, (half_height_r1 - d) * -j_cap, std_dir, SC.Manifold(r_c1 + d, -j_cap, (r_c1 - (r_c1 + d) / 2) * -j_cap)),
            (c1, r1, (half_height_r1 - d) * j_cap, std_dir, SC.Manifold(r_c1 + d, j_cap, (r_c1 - (r_c1 + d) / 2) * j_cap)),
            (c1, r1, (half_height_r1 + d) * j_cap, std_dir, SC.Manifold(r_c1 - d, j_cap, (r_c1 - (r_c1 - d) / 2) * j_cap)),
            (c1, r1, (half_height_r1 + r_c1 - d) * j_cap, std_dir, SC.Manifold(d, j_cap, (r_c1 - d / 2) * j_cap)),

            (c1, r1, top_right_r1 + (r_c1 - d) * unit_45, std_dir, SC.Manifold(d, unit_45, (r_c1 - d / 2) * unit_45)),

            # reverse check with rotated_dir
            (c1, r1, (r_c1 - d) * unit_45 + SC.rotate(top_right_r1, rotated_dir), rotated_dir, SC.Manifold(d, unit_45, (r_c1 - d / 2) * unit_45)),
            ]

            test_manifold_list(manifold_list)
        end

        Test.@testset "Rect2D vs. Rect2D" begin
            manifold_list_no_dir = [
            # std_dir
            (r2, r1, (half_height_r2 + half_height_r1 - d) * -j_cap, SC.Manifold(d, -j_cap, (half_height_r2 - d/2) * -j_cap)),
            (r2, r1, (half_height_r2 + d) * -j_cap, SC.Manifold(half_height_r1 - d, -j_cap, (half_height_r2 - (half_height_r1 - d)/2) * -j_cap)),
            (r2, r1, (half_height_r2 - d) * -j_cap, SC.Manifold(half_height_r1 + d, -j_cap, (half_height_r2 - (half_height_r1 + d)/2) * -j_cap)),
            (r2, r1, d * -j_cap, SC.Manifold(half_height_r2 + half_height_r1 - d, -j_cap, d * -j_cap)),
            (r2, r1, d * j_cap, SC.Manifold(half_height_r2 + half_height_r1 - d, j_cap, d * j_cap)),
            (r2, r1, (half_height_r2 - d) * j_cap, SC.Manifold(half_height_r1 + d, j_cap, (half_height_r2 - (half_height_r1 + d)/2) * j_cap)),
            (r2, r1, (half_height_r2 + d) * j_cap, SC.Manifold(half_height_r1 - d, j_cap, (half_height_r2 - (half_height_r1 - d)/2) * j_cap)),
            (r2, r1, (half_height_r2 + half_height_r1 - d) * j_cap, SC.Manifold(d, j_cap, (half_height_r2 - d/2) * j_cap)),

            (r2, r1, (half_width_r2 + half_width_r1 - d) * -i_cap, SC.Manifold(d, -i_cap, (half_width_r2 - d/2) * -i_cap)),
            (r2, r1, (half_width_r2 + d) * -i_cap, SC.Manifold(half_width_r1 - d, -i_cap, (half_width_r2 - (half_width_r1 - d)/2) * -i_cap)),
            (r2, r1, (half_width_r2 - d) * -i_cap, SC.Manifold(half_width_r1 + d, -i_cap, (half_width_r2 - (half_width_r1 + d)/2) * -i_cap)),
            (r2, r1, (half_width_r2 - d) * i_cap, SC.Manifold(half_width_r1 + d, i_cap, (half_width_r2 - (half_width_r1 + d)/2) * i_cap)),
            (r2, r1, (half_width_r2 + d) * i_cap, SC.Manifold(half_width_r1 - d, i_cap, (half_width_r2 - (half_width_r1 - d)/2) * i_cap)),
            (r2, r1, (half_width_r2 + half_width_r1 - d) * i_cap, SC.Manifold(d, i_cap, (half_width_r2 - d/2) * i_cap)),

            (r2, r1, top_right_r2 .- d, SC.Manifold(half_height_r1 + d, j_cap, top_right_r2 + (half_width_r1 + d)/2 * -i_cap + (half_height_r1 + d)/2 * -j_cap)),
            (r2, r1, top_right_r2, SC.Manifold(half_height_r1, j_cap, top_right_r2 + (half_width_r1/2) * -i_cap + (half_height_r1/2) * -j_cap)),
            (r2, r1, top_right_r2 .+ d, SC.Manifold(half_height_r1 - d, j_cap, top_right_r2 + (half_width_r1 - d)/2 * -i_cap + (half_height_r1 - d)/2 * -j_cap)),
            ]

            test_manifold_list_no_dir(manifold_list_no_dir)

            manifold_list = [
            # std_dir
            (r2, r1, (half_height_r2 + half_height_r1 - d) * -j_cap, std_dir, SC.Manifold(d, -j_cap, (half_height_r2 - d/2) * -j_cap)),
            (r2, r1, (half_height_r2 + d) * -j_cap, std_dir, SC.Manifold(half_height_r1 - d, -j_cap, (half_height_r2 - (half_height_r1 - d)/2) * -j_cap)),
            (r2, r1, (half_height_r2 - d) * -j_cap, std_dir, SC.Manifold(half_height_r1 + d, -j_cap, (half_height_r2 - (half_height_r1 + d)/2) * -j_cap)),
            (r2, r1, d * -j_cap, std_dir, SC.Manifold(half_height_r2 + half_height_r1 - d, -j_cap, d * -j_cap)),
            (r2, r1, d * j_cap, std_dir, SC.Manifold(half_height_r2 + half_height_r1 - d, j_cap, d * j_cap)),
            (r2, r1, (half_height_r2 - d) * j_cap, std_dir, SC.Manifold(half_height_r1 + d, j_cap, (half_height_r2 - (half_height_r1 + d)/2) * j_cap)),
            (r2, r1, (half_height_r2 + d) * j_cap, std_dir, SC.Manifold(half_height_r1 - d, j_cap, (half_height_r2 - (half_height_r1 - d)/2) * j_cap)),
            (r2, r1, (half_height_r2 + half_height_r1 - d) * j_cap, std_dir, SC.Manifold(d, j_cap, (half_height_r2 - d/2) * j_cap)),

            (r2, r1, (half_width_r2 + half_width_r1 - d) * -i_cap, std_dir, SC.Manifold(d, -i_cap, (half_width_r2 - d/2) * -i_cap)),
            (r2, r1, (half_width_r2 + d) * -i_cap, std_dir, SC.Manifold(half_width_r1 - d, -i_cap, (half_width_r2 - (half_width_r1 - d)/2) * -i_cap)),
            (r2, r1, (half_width_r2 - d) * -i_cap, std_dir, SC.Manifold(half_width_r1 + d, -i_cap, (half_width_r2 - (half_width_r1 + d)/2) * -i_cap)),
            (r2, r1, (half_width_r2 - d) * i_cap, std_dir, SC.Manifold(half_width_r1 + d, i_cap, (half_width_r2 - (half_width_r1 + d)/2) * i_cap)),
            (r2, r1, (half_width_r2 + d) * i_cap, std_dir, SC.Manifold(half_width_r1 - d, i_cap, (half_width_r2 - (half_width_r1 - d)/2) * i_cap)),
            (r2, r1, (half_width_r2 + half_width_r1 - d) * i_cap, std_dir, SC.Manifold(d, i_cap, (half_width_r2 - d/2) * i_cap)),

            (r2, r1, top_right_r2 .- d, std_dir, SC.Manifold(half_height_r1 + d, j_cap, top_right_r2 + (half_width_r1 + d)/2 * -i_cap + (half_height_r1 + d)/2 * -j_cap)),
            (r2, r1, top_right_r2, std_dir, SC.Manifold(half_height_r1, j_cap, top_right_r2 + (half_width_r1/2) * -i_cap + (half_height_r1/2) * -j_cap)),
            (r2, r1, top_right_r2 .+ d, std_dir, SC.Manifold(half_height_r1 - d, j_cap, top_right_r2 + (half_width_r1 - d)/2 * -i_cap + (half_height_r1 - d)/2 * -j_cap)),

            # rotated_dir
            (r2, r1, (half_height_r2 - d) * -j_cap - SC.rotate(top_right_r1, rotated_dir), rotated_dir, SC.Manifold(d, -j_cap, ((zero(T) - d / tan(theta) + d * tan(theta)) / 3) * i_cap + ((-half_height_r2 + d - half_height_r2 - half_height_r2) / 3) * j_cap)),
            (r2, r1, (half_height_r2 - d) * j_cap + SC.rotate(top_right_r1, rotated_dir), rotated_dir, SC.Manifold(d, j_cap, ((zero(T) + d / tan(theta) - d * tan(theta)) / 3) * i_cap + ((half_height_r2 - d + half_height_r2 + half_height_r2) / 3) * j_cap)),

            (r2, r1, (half_width_r2 + LA.norm(top_right_r1) * cos(-theta_r1 + theta) - d) * -i_cap, rotated_dir, SC.Manifold(d, -i_cap, ((-half_width_r2 + d - half_width_r2 - half_width_r2) / 3) * i_cap + (LA.norm(top_right_r1) * sin(-theta_r1 + theta) + (d / tan(theta) - d * tan(theta)) / 3) * j_cap)),
            (r2, r1, (half_width_r2 + LA.norm(top_right_r1) * cos(-theta_r1 + theta) - d) * i_cap, rotated_dir, SC.Manifold(d, i_cap, ((half_width_r2 - d + half_width_r2 + half_width_r2) / 3) * i_cap + (LA.norm(top_right_r1) * sin(convert(T, pi - theta_r1) + theta) + (d * tan(theta) - d / tan(theta)) / 3) * j_cap)),
            ]

            test_manifold_list(manifold_list)
        end
    end
end
