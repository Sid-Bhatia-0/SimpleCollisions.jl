# PhysicsPrimitives2D

This package provides **lightweight** and **efficient** primitives for 2D physics. This package is inspired by [ImpulseEngine](https://github.com/RandyGaul/ImpulseEngine) and [box2d](https://github.com/erincatto/box2d).

# Design
1. [Geometry](#geometry)
    1. [Axes](#axes)
    1. [Shapes](#shapes)
        1. [StdPoint](#stdpoint)
        1. [StdLine](#stdline)
        1. [StdCircle](#stdcircle)
        1. [StdRect](#stdrect)
1. [Collisions](#collisions)
    1. [Collision Detection](#collision-detection)
    1. [Collision Manifold](#collision-manifold)
1. [Physics](#physics)

## Geometry

### Axes

A position on the 2D coordinate plane (with respect to the world origin) is represented using an instance of `GeometryBasics.Vec{2, T}`.

Orientation on the 2D coordinate plane (with respect to the world x-axis) can be represented compactly by a single scalar denoting the counterclockwise angle `theta` (in radians) with respect to the x-axes. Equivalently, it can be represented by a 2D unit vector with components `cos(theta)` and `sin(theta)`. We represent orientation with the following struct:

```
struct Axes{T}
    x_cap::GB.Vec{2, T}
    y_cap::GB.Vec{2, T}
end
```

Here `x_cap` (representing the local x-axis) and `y_cap` (representing the local y-axis) are supposed to be orthogonal unit vectors. This representation allows for caching of `cos(theta)` and `sin(theta)`, thereby avoiding unnecessary compution. We could have used just stored `x_cap`. Additionally storing `y_cap` is just a matter of taste and has minimal overhead.

### Shapes

This package provides the following shapes: `StdPoint`, `StdLine`, `StdCircle`, and `StdRect`.
the `Std` prefix stands for the word standard. Here, a standard shape refers to a shape whose geometric center is at the world origin and whose axes of orientation (determined by symmetry) are aligned with the world coordinate axes.

A standard shape can be augmented with a position and (optionally) an axes to represent that shape at an arbitrary location and orientation with respect to the world frame of reference. This decoupling of the shape, position, and orientation is useful for collision detection and collision manifold generation as it allows us to exploit the symmetry offered in a use case.

#### StdPoint

```
struct StdPoint{T} <: AbstractStdShape{T} end
```

`StdPoint` is a point placed at the origin. It doesn't require any fields.

<img src="https://github.com/Sid-Bhatia-0/PhysicsPrimitives2D.jl/raw/master/docs/assets/img/StdPoint.svg" width="360px">

#### StdLine

```
struct StdLine{T} <: AbstractStdShape{T}
    half_length::T
end
```

`StdLine` is a line segment aligned with the world x-axes centered at the origin. It requires only one field - a `half_length`.

<img src="https://github.com/Sid-Bhatia-0/PhysicsPrimitives2D.jl/raw/master/docs/assets/img/StdLine.svg" width="360px">

#### StdCircle

```
struct StdCircle{T} <: AbstractStdShape{T}
    radius::T
end
```

`StdCirle` is a circle centered at the origin. It requires only one field - a `radius`.

<img src="https://github.com/Sid-Bhatia-0/PhysicsPrimitives2D.jl/raw/master/docs/assets/img/StdCircle.svg" width="360px">

#### StdRect

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

This package offers the `collision_detection` function to detect collisions between pairs of standard shapes with arbitrary relative location and orientation.

Collision detection in this package relies on relative position and orientation rather than absolute ones. This is to facilitate exploitation of symmetry for some collision detection computations that can be made more efficient. For examples, in the case of `StdRect` vs. `StdRect`, it is much easier to detect a collision and generate a collision manifold if we know that both shapes are axes aligned with respect to the world coordinates.

### Collision Manifold

A collision manifold contains information about how to resolve a collision once it has been detected. The following is the struct used to represent the collision manifold:

```
struct Manifold{T}
    penetration::T
    axes::Axes{T}
    contact::GB.Vec{2, T}
end
```

Here the `axes` field provides the collision tangent and normal. The `x_cap` field of this `axes` corresponds to the collision tangent and the `y_cap` field corresponds to the collision normal. A `Manifold` object is calculated relative to a standard shape.

## Physics

This package provides the following function to calculate the linear impulse for two collising rigid bodies:
1. `get_normal_impulse` function calcuates the normal impulse for the following three types of collisions:
    1. Dynamic body vs. Dynamic body
    1. Dynamic body vs. Kinetic body
    1. Dynamic body vs. Static body
1. `get_tangential_impulse` function calcuates the tangential impulse (caused by friction)
