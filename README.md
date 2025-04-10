# DiscretePhasePortrait
The DiscretePhasePortrait package provides one function:

	phaseportrait(F, G, limits; <keyword arguments>)

Plot a phase portrait for a discrete 2D system.

A phase portrait of the discrete-time system

```math
	x_{t+1} =  F(x_t, y_t),     y_{t+1} = G(x_t, y_t)
```

is drawn showing isoclines at particular directions and their preimages for the rectangular 
region of the plane defined by `limits=[xmin xmax; ymin ymax]`.  Labels and a legend are used 
to indicate which regions map to others and the direction of motion in each region.

# Keyword Arguments

  - `fixpt`: Approximation `(x0,y0)` of a fixed point of interest.
  - `Jac`: Function of `x` and `y` returning the Jacobian matrix of the system at a point `(x,y)`.
  - `directions=[0,π/2]`: Directions of the isoclines to use if `fixpt` or `Jac` are not given.

> [!NOTE]
> If both `fixpt` and `Jac` are specified, the actual fixed point is determined using
> Newton's method and `Jac` starting from the input `fixpt`, then the eigenvalues and
> eigenvectors of the Jacobian at the actual fixed point are calculated and the
> directions of the isoclines are set to the eigenvector directions.  The return
> value is `(actualfixedpoint,eigenvalues,eigenvectors,evdirections)`.  If either
> `fixpt` or `Jac` are *not* specified then the directions of the isoclines are
> determined by the values in `directions`, and the return value is `nothing`.

  - `n::Int=1`: Number of iterations of the function f=[F,G] applied.
  - `labelpoints`: User defined points at which to apply region labeling.  

    This may be an array or tuple of `(x,y)` points.  If `labelpoints` is not specified then 
    by default, four points on the main plot near the four corners are labelled.

  - `showdet0::Bool=true`: Plot the det(Jac)=0 curve and its image. `Jac` must also be given.
  - `detlimits=limits`: Limits of the form `[dxmin dxmax; dymin dymax]` to use for the image of the 
    det(Jac)=0 curve.
    The points within `detlimits` satisfying det(Jac)=0 are computed and the image of all such
    points are plotted if they lie within `limits`. 

> [!TIP]
> Often, to see more of the image curve it is 
> necessary to make `detlimits` a larger rectangle than `limits`.

  - `showrange::Bool=true`: Gray shade the image of all points within `limits`.
  - `resolution::Int=500`: Grid resolution.
  - `plotsize=(800,600)`: Set the size of the plot in pixels.


[![Build Status](https://github.com/allanwillms/DiscretePhasePortrait.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/allanwillms/DiscretePhasePortrait.jl/actions/workflows/CI.yml?query=branch%3Amain)
