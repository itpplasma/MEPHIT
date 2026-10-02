# Maxwell exterior boundary

MEPHIT solves the magnetic response on the plasma mesh and a surrounding
vacuum annulus. `config%maxwell_outer_scale` sets the semiaxes of the annulus's
outer ellipse relative to the plasma bounding box. Its default, `2`, retains
the existing coordinates, point count and ordering. Increasing the value
allows a numerical study of the exterior boundary while retaining the plasma
mesh and forcing.

The unchanged FreeFEM boundary condition sets tangential vector-potential
edge degrees of freedom to zero. It therefore imposes zero normal magnetic
response at this finite boundary. Changing its distance is a numerical control;
it does not introduce a free-space boundary condition or establish accuracy.
Exterior mesh convergence must be checked alongside boundary distance.

The scale must be finite and exceed one. Before triangulation, the generated
outer polygon must lie at positive cylindrical radius and strictly contain
every plasma boundary vertex. Containment is checked against the polygon's
chords; checking only an analytic ellipse would miss some small-scale failures.
Invalid requests are rejected without clamping.

The effective value is stored under `/mesh/maxwell_outer_scale`. Cached mesh
replay checks this immutable geometry value against the requested configuration
before reading mesh arrays or initializing FreeFEM. Legacy meshes without the
dataset have effective scale `2`. A mismatch requires remeshing and rebuilding
the preconditioner; overwriting `/config` does not change the stored geometry.
