import .Makie

get_angle_reference(angle::T) where {T} = GB.Point2{T}(cos(angle), sin(angle))

function get_angle_reference(shape::GB.Circle, angle)
    center = get_center(shape)
    return [center, center .+ shape.r .* get_angle_reference(angle)]
end

get_angle_reference(body::RigidBody) = get_angle_reference(get_shape(body), get_angle(body))

function init_screen(world_node::Makie.Observable{<:World}; resolution = (720, 720), xlims = (0, 10), ylims = (0, 10))
    scene = Makie.Scene(resolution = resolution)

    shapes_node = Makie.lift(world_node) do world
        bodies = get_bodies(world)
        shapes = get_shape.(bodies)
        return shapes
    end

    Makie.poly!(scene, shapes_node, color = :lightgray, linestyle = :solid, strokewidth = 2.0, strokecolor = :black)

    angle_references_node = Makie.lift(world_node) do world
        bodies = get_bodies(world)
        angle_references = get_angle_reference.(bodies)
        return angle_references
    end

    Makie.poly!(scene, angle_references_node, color = :lightgray, linestyle = :solid, strokewidth = 2.0, strokecolor = :black)

    Makie.xlims!(scene, xlims)
    Makie.ylims!(scene, ylims)

    return scene
end

function render!(world::World; num_iter = 100, dt = 0.01, frame_rate = 24)
    world_node = Makie.Node(world)
    scene = init_screen(world_node)
    Makie.display(scene)
    for i in 1:num_iter
        world_node[] = step!(world_node[], dt)
        sleep(1 / frame_rate)
    end
    return scene
end
