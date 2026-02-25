# MultiStreamVlasovPoisson.jl

Documentation for MultiStreamVlasovPoisson.jl


## Installation

```@example main
using MultiStreamVlasovPoisson
using Plots

eps = 1.0
nx = 200
k = 0.5
xmin, xmax = 0.0, 2π / k
vmin, vmax = -6.0, 6.0
ng = 200
mesh = UniformMesh(xmin, xmax, nx, vmin, vmax, ng)

rho, u, rho_tot = compute_initial_condition(mesh, k)

poisson = NonLinearPoissonSolver(eps, nx)

phi = -log.(rho_tot)
plot(mesh.x, rho_tot)
```

## Functions

```@autodocs
Modules = [MultiStreamVlasovPoisson]
Order   = [:type, :function]
```
