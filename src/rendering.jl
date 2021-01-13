import .Makie

function init_screen(world_node::Makie.Observable{<:World}; resolution = (720, 720))
    scene = Makie.Scene(resolution = resolution)
    Makie.poly!(scene, [body.shape for body in world_node[].bodies])
    return scene
end

function render(world::World)
    world_node = Makie.Node(world)
    scene = init_screen(world_node)
    Makie.display(scene)
end
