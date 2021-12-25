# SimpleCollisions

This package provides **fast** and **lightweight** "primitives" for 2D physics. Roughly speaking, it offers a collection of methods and structs needed to detect and resolve collisions between simple 2D shapes like circles and rectangles, along with some basic methods to calculate impulses when two rigid bodies collide with each other. The priority is high performance for the simple cases instead of adding sophisticated feature.

In order to cite this package, please refer to the file `CITATION.bib`. Starring the repository on GitHub is also appreciated.

# Design
1. [What This Package Is NOT About](#what-this-package-is-not-about)
1. [Position and Orientation](#position-and-orientation)
1. [Shapes](#shapes)
    1. [StdPoint](#stdpoint)
    1. [StdLine](#stdline)
    1. [StdCircle](#stdcircle)
    1. [StdRect](#stdrect)
1. [Collisions](#collisions)
    1. [Collision Detection](#collision-detection)
    1. [Collision Manifold](#collision-manifold)
1. [Physics](#physics)
    1. [Impulse](#impulse)

## What This Package Is NOT About

This package is not a physics *engine* in itself, even though it is inspired by physics engines [ImpulseEngine](https://github.com/RandyGaul/ImpulseEngine) and [box2d](https://github.com/erincatto/box2d). The reason for only providing primitives instead of a full-blown physics engine is to keep this package lightweight, make minimal or no assumptions about a game, and allow for context-specific optimizations in a game (that are often hard to generalize in advance).

This package does not provide a generalized game loop. In order to understand why this is so, consider the example of 2D collision detection. In *your* game loop, it may or may not be terribly inefficient to check for exact pair-wise collisions (narrow phase) for all possible pairs of rigid bodies at every iteration. And there are several possible heuristics to optimize this further (like broad phase algorithms). *But*, once you know the specific game you are trying to build, you get full control of the collision detection process and can choose to check only those pairs that are make sense in the context of your game. This not only allows for a fine degree of optimization, but could also be a much simpler solution overall as compared to the generalized broad-phase heuristics that may still not be good enough for your case.

## Position and Orientation

A position on the 2D coordinate plane (with respect to the world origin) is represented using an instance of `SC.Vector2D{T}`. An orientation on the 2D coordinate plane (with respect to the world coordinate axes) can be represented compactly by a single scalar denoting the counterclockwise angle `theta` (in radians) with respect to the world x-axes. Equivalently, a 2D unit vector `SC.Vector2D(cos(theta), sin(theta))` could also be used to represent a direction. This package uses the latter representation of orientation.

Computations involving orientation often require `cos(theta)` and `sin(theta)`. Overall, it is cheaper to cache `cos(theta)`s and `sin(theta)`s than to repeatedly compute the trigonmetric functions from the angle `theta` each time.

## Shapes

This package provides the following shapes: `StdPoint`, `StdLine`, `StdCircle`, and `StdRect`. The prefix `Std` stands for standard. Here, a standard shape refers to a shape whose geometric center is placed at the world origin and whose orientation (as determined by some form of symmetry) is aligned with the positive world x-axis. We use `Std` shapes to decouple the shape of a body from its position and orientation with respect to the world frame of reference. A `Std` shape can be passed along with a position and (optionally) an orientation to represent a body of that shape at an arbitrary location and orientation with respect to the world frame of reference.

### StdPoint

```
struct StdPoint{T} <: AbstractStdShape{T} end
```

`StdPoint` is a point placed at the origin. It doesn't require any fields.

<img src="https://github.com/Sid-Bhatia-0/SimpleCollisions.jl/raw/master/docs/assets/img/StdPoint.svg" width="360px">

### StdLine

```
struct StdLine{T} <: AbstractStdShape{T}
    half_length::T
end
```

`StdLine` is a line segment centered at the origin and aligned with the world x-axes. It requires only one field - a `half_length`.

<img src="https://github.com/Sid-Bhatia-0/SimpleCollisions.jl/raw/master/docs/assets/img/StdLine.svg" width="360px">

### StdCircle

```
struct StdCircle{T} <: AbstractStdShape{T}
    radius::T
end
```

`StdCirle` is a circle centered at the origin. It requires only one field - a `radius`.

<img src="https://github.com/Sid-Bhatia-0/SimpleCollisions.jl/raw/master/docs/assets/img/StdCircle.svg" width="360px">

### StdRect

```
struct StdRect{T} <: AbstractStdShape{T}
    half_width::T
    half_height::T
end
```

`StdRect` is a rectangle centered at the origin with its edges parallel to the world coordinate axes (width edges parallel to world x-axis and height edges parallel to world y-axis). It requires two fields - a `half_width` and a `half_height`.

<img src="https://github.com/Sid-Bhatia-0/SimpleCollisions.jl/raw/master/docs/assets/img/StdRect.svg" width="360px">

## Collisions

### Collision Detection

This package offers the `is_colliding` function to detect collisions between pairs of bodies with `Std` shapes at arbitrary relative positions and orientations. If the orientation of both the bodies are aligned with world x-axis (for example, collision detection between two axes-aligned bounding boxes), that is, the bodies are just shifted versions of `Std` shapes, or even if their orientations only mutually align, then do not pass the relative orientation argument to the `is_colliding` function and you will dispatch to faster method.

The position and orientation arguments (whenever present) in the `is_colliding` method definitinos refer to the relative position and orientation of the second body (second argument) in the frame of reference of the first body (first argument). The decoupling of the shapes of the bodies shapes from their positions and orientations allows us to exploit symmetry and speed up the computation for collision detection and manifold generation for certain common use cases (for example, in the case of collision of two axes-aligned bounding boxes, where you would not pass the relative orientation argument to `is_colliding` function). We use relative positions and orientations instead of absolute ones in order to make the geometry calculations easier to understand. I don't think frame of reference conversions required for this have any significant computational overhead. Calculating relative positions is just subtracting two vectors, and calculating relative orientation is pretty cheap too because we already cache `sin(theta)` and `cos(theta)` in the unit vectors.

### Collision Manifold

A collision manifold contains information about how to resolve a collision once it has been detected. The following is the struct used to represent a collision manifold:

```
struct Manifold{T}
    penetration::T
    normal::Vector2D{T}
    contact::Vector2D{T}
end
```

Here the `normal` field corresponds to the collision normal. A `Manifold` object is calculated with respect to the frame of reference of the first body (first argument). **It is important to note that manifold generation for a collision assumes (wherever needed) that the two objects are indeed colliding.**

## Physics

This package does not provide a concrete type representing a generalized rigid body. Some rigid bodies may be static (may not require any mass or velocity related fields), some may be kinetic (may not require mass related fields), some may be dynamic (may require a lot of fields), some may not rotate at all (may not require any orientation related fields), some may be frictionless (may not require static or kinematic friction coefficients) etc... Basically, there is a lot of variety out there, and trying to make an optimized decision upfront is very difficult, if not impossible. Storing all such potential fields into one rigid body struct would often lead to wasteful memory and computations.

This package provides the following functions to calculate the linear impulse of two colliding rigid bodies:
1. `get_normal_impulse`: calcuates the normal impulse for a collisions
1. `get_tangential_impulse`: calcuates the tangential impulse (caused by friction) for a collision
