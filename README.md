# PhysicsPrimitives2D

This package provides **lightweight** and **efficient** primitives for 2D physics. This package is inspired by [ImpulseEngine](https://github.com/RandyGaul/ImpulseEngine) and [box2d](https://github.com/erincatto/box2d).

# Design
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

## Position and Orientation

A position on the 2D coordinate plane (with respect to the world origin) is represented using an instance of `StaticArrays.SVector{2, T}`. An orientation on the 2D coordinate plane (with respect to the world coordinate axes) can be represented compactly by a single scalar denoting the counterclockwise angle `theta` (in radians) with respect to the world x-axes. Equivalently, a 2D unit vector `StaticArrays.SVector(cos(theta), sin(theta))` could also be used to represent a direction. Finally, one could also store the 2D [rotation matrix](https://en.wikipedia.org/wiki/Rotation_matrix) corresponding to the angle `theta`. This package uses the rotation matrix to represent orientation. We store the two orthogonal unit vectors of the rotation matrix inside the `Axes` struct:

```
struct Axes{T}
    x_cap::SA.SVector{2, T}
    y_cap::SA.SVector{2, T}
end
```

Computations involving orientation often require `cos(theta)` and `sin(theta)`. Overall, it is cheaper to cache `cos(theta)`s and `sin(theta)`s in the `Axes` struct than to repeatedly compute the trigonmetric functions from the angle `theta` each time.

## Shapes

This package provides the following shapes: `StdPoint`, `StdLine`, `StdCircle`, and `StdRect`.
the `Std` prefix stands for the word standard. Here, a standard shape refers to a shape whose geometric center is at the world origin and whose axes of orientation (determined by symmetry) are aligned with the world coordinate axes.

A standard shape can be augmented with a position and (optionally) an axes to represent that shape at an arbitrary location and orientation with respect to the world frame of reference. This decoupling of shape, position, and orientation of a body allows us to exploit symmetry and speed up the computations for collision detection and collision manifold generation for certain common use cases (for example, collision of two axes-aligned bounding boxes).

### StdPoint

```
struct StdPoint{T} <: AbstractStdShape{T} end
```

`StdPoint` is a point placed at the origin. It doesn't require any fields.

<img src="https://github.com/Sid-Bhatia-0/PhysicsPrimitives2D.jl/raw/master/docs/assets/img/StdPoint.svg" width="360px">

### StdLine

```
struct StdLine{T} <: AbstractStdShape{T}
    half_length::T
end
```

`StdLine` is a line segment aligned with the world x-axes centered at the origin. It requires only one field - a `half_length`.

<img src="https://github.com/Sid-Bhatia-0/PhysicsPrimitives2D.jl/raw/master/docs/assets/img/StdLine.svg" width="360px">

### StdCircle

```
struct StdCircle{T} <: AbstractStdShape{T}
    radius::T
end
```

`StdCirle` is a circle centered at the origin. It requires only one field - a `radius`.

<img src="https://github.com/Sid-Bhatia-0/PhysicsPrimitives2D.jl/raw/master/docs/assets/img/StdCircle.svg" width="360px">

### StdRect

```
struct StdRect{T} <: AbstractStdShape{T}
    half_width::T
    half_height::T
end
```

`StdRect` is a rectangle centered at the origin with its edges parallel to the world coordinate axes. It requires two fields - a `half_width` and a `half_height`.

<img src="https://github.com/Sid-Bhatia-0/PhysicsPrimitives2D.jl/raw/master/docs/assets/img/StdRect.svg" width="360px">

## Collisions

### Collision Detection

This package offers the `is_colliding` function to detect collisions between pairs of standard shapes with arbitrary relative positions and orientations. If the shapes are axes aligned (for example, collision detection between two axes-aligned bounding boxes), that is, they are just shifted versions of `Std` shapes, then do not pass the relative axes argument to the `is_colliding` function and you will dispatch to faster method. Exploitation of symmetry is the reason why decouples shape, position, and orientation.

The position and axes arguments (whenever present) in the `is_colliding` function are the relative position and orientation of the second shape (second argument) with respect to the frame of reference of the first shape (first argument).

### Collision Manifold

A collision manifold contains information about how to resolve a collision once it has been detected. The following is the struct used to represent a collision manifold:

```
struct Manifold{T}
    penetration::T
    axes::Axes{T}
    contact::SA.SVector{2, T}
end
```

Here the `axes` field contains the collision tangent and normal. The `x_cap` field of this `axes` corresponds to the collision tangent and the `y_cap` field corresponds to the collision normal. A `Manifold` object is calculated with respect to the frame of reference of the first object (first argument). **It is important to note that collision manifold generation assumes (wherever required) that the two objects are indeed colliding.**

## Physics

This package provides the following functions to calculate the linear impulse of two colliding rigid bodies:
1. `get_normal_impulse`: calcuates the normal impulse for a collisions
1. `get_tangential_impulse`: calcuates the tangential impulse (caused by friction) for a collision
