# Closed equilibrium boundary

For a smooth equilibrium with nested closed flux surfaces, enable an explicit
plasma boundary in the `scalars` namelist:

```fortran
config%closed_boundary_psi = 0.9999
```

The value is the fraction of the EQDSK axis-to-boundary poloidal flux and must
be in `(0,1]`. A negative value retains the existing diverted-equilibrium
boundary search. The selected contour supplies the plasma geometry; changing
the EQDSK rectangle does not deliberately change that contour. The vacuum
Maxwell boundary remains a separate numerical boundary.

This mode requires a single outward flux crossing on the magnetic-axis
midplane and nested contours that fit inside the EQDSK rectangle. The root
search brackets a sampled outward crossing; it cannot resolve arbitrary pairs
of roots between samples. It is unsuitable for multiple magnetic axes,
non-nested surfaces, or a separatrix with an X point. Coordinate rescaling via
`circ_mesh_scale` is incompatible with this option.

The surface tracer rejects ODE evaluations outside the rectangle, including
trial stages, and checks contour closure and the
requested flux. Analytic tests cover both flux signs, flux offsets, vertical
axis shifts and rectangle sizes, and compare the native safety factor, area,
perimeter and boundary coordinates with exact circular contours.

A contour very close to a rectangle face can be rejected when an adaptive ODE
trial stage leaves the supplied field domain. Supply adequate EQDSK padding;
the closed contour remains defined by the requested flux.
