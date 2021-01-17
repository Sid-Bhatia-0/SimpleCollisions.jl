import .Makie

function init_screen(world_node::Makie.Observable{<:World}; resolution = (720, 720), xlims = (0, 10), ylims = (0, 10))
    scene = Makie.Scene(resolution = resolution)
    shapes_node = Makie.lift(world_node) do world
        shapes = get_shape.(get_bodies(world))
    end
    Makie.poly!(shapes_node)
    Makie.xlims!(scene, xlims)
    Makie.ylims!(scene, ylims)
    return scene
end

function render!(world::World; num_iter = 100, dt = 0.01, frame_rate = 24)
    world_node = Makie.Node(world)
    scene = init_screen(world_node)
    Makie.display(scene)
    for i in 1:num_iter
        world_node[] = update!(world_node[], dt)
        sleep(1 / frame_rate)
    end
    return scene
end
