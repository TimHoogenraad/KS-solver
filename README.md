# Kuramoto-Sivashinsky solver

A small MATLAB solver for the one-dimensional Kuramoto-Sivashinsky (KS)
equation. It supports periodic boundaries and zero-valued (homogeneous
Dirichlet) boundaries, and returns the complete solution history for plotting
or analysis.

The KS equation is a standard model for nonlinear pattern formation and
spatiotemporal chaos. In the convention used here, it is

```text
u_t + u u_x + v2 u_xx + v4 u_xxxx = 0.
```

The second-derivative term destabilizes long waves, the fourth-derivative term
damps short waves, and the nonlinear term transfers energy between scales. The
default values are `v2 = 1` and `v4 = 1`.

## Requirements

- MATLAB
- No additional MATLAB toolboxes

Clone or download this repository, start MATLAB, and change to the repository
directory. Add the solver to the MATLAB path:

```matlab
addpath('src')
```

## Quick start

Create a configuration structure and pass it to `solveKS`:

```matlab
config.boundary = 'periodic';
config.L = 64;       % Domain length
config.N = 128;      % Number of grid points
config.dt = 0.01;    % Time-step size
config.steps = 100;  % Number of time steps

config.initial = @(x) 0.5*sin(2*pi*x/config.L) ...
    .* (1 + 0.3*sin(4*pi*x/config.L));

[t,x,u] = solveKS(config);

imagesc(t,x,u)
axis xy
xlabel('Time')
ylabel('x')
colorbar
```

Here, `u(i,j)` is the solution at position `x(i)` and time `t(j)`. The first
column, `u(:,1)`, contains the initial condition.

## Boundary conditions

Set `config.boundary` to one of the following values:

| Value | Meaning | Grid size returned in `x` |
| --- | --- | --- |
| `'periodic'` | The left and right sides of the domain connect. | `N` |
| `'dirichlet'` | The solution is fixed to zero at `x = 0` and `x = L`. | `N + 2` |

For a periodic run, `N` must be even. The grid contains
`L/N, 2L/N, ..., L`; the points at `0` and `L` represent the same periodic
location, so only `L` is stored.

For a Dirichlet run, `N` is the number of interior points. The returned grid
also includes the two boundary points, giving `N + 2` values in total. The
solver sets both endpoint values to zero, including in the initial condition.

To switch the quick-start example to zero-valued boundaries, change only:

```matlab
config.boundary = 'dirichlet';
```

## Configuration reference

| Field | Required | Description |
| --- | --- | --- |
| `boundary` | Yes | Either `'periodic'` or `'dirichlet'`. |
| `L` | Yes | Positive domain length. |
| `N` | Yes | Periodic grid points, or Dirichlet interior grid points. Must be an integer of at least 2; periodic runs require an even value. |
| `dt` | Yes | Positive time-step size. |
| `steps` | Yes | Nonnegative integer number of time steps. |
| `initial` | Yes | Function handle evaluated on `x`, or a numeric vector with one value for each returned grid point. |
| `v2` | No | Coefficient of `u_xx`. Defaults to `1`. |
| `v4` | No | Coefficient of `u_xxxx`. Defaults to `1`. |
| `Cs` | No | Nonnegative subgrid-scale model coefficient. If omitted, the model is disabled. |

When `initial` is a function handle, it must accept the grid vector and return
one value per grid point. For example:

```matlab
config.initial = @(x) sin(2*pi*x/config.L);
```

A numeric initial condition is also accepted:

```matlab
x0 = config.L*(1:config.N)'/config.N;  % Periodic grid
config.initial = sin(2*pi*x0/config.L);
```

For a Dirichlet run, a numeric initial condition must contain `N + 2` values,
including the two endpoints. Any nonzero endpoint values are replaced by zero.

## Outputs

```matlab
[t,x,u] = solveKS(config);
```

| Output | Shape | Description |
| --- | --- | --- |
| `t` | `1 x (steps + 1)` | Times from `0` through `steps*dt`. |
| `x` | `N x 1` or `(N + 2) x 1` | Spatial grid. |
| `u` | `numel(x) x (steps + 1)` | Complete solution history. |

Because `solveKS` stores every state, memory use grows with both `N` and
`steps`. For example, doubling either value approximately doubles the memory
required for `u`.

## Included example

The example script runs the same initial condition with both boundary types
from `t = 0` through `t = 500`. This includes the initial transient and a
window of developed chaotic dynamics. It displays the solutions as stacked
grayscale space-time surfaces using the lighting style of the original KS
visualization code:

```matlab
run('examples/compare_boundaries.m')
```

The script locates `src` automatically, even when it is launched from another
MATLAB working directory using its absolute path.

The same original-style visualization can be used for any solution returned
by `solveKS`:

```matlab
plotKSSpacetime(t,x,u,'DomainLength',config.L,'Title','Periodic')
```

`plotKSSpacetime` creates a figure by default. Use the optional `Parent` value
to draw into an existing axes, or `ColorLimit` to give several plots the same
symmetric color scale.

## Numerical method

The solver uses a Fourier spectral discretization in space and the fourth-order
exponential time-differencing Runge-Kutta method (ETDRK4) in time.

- Periodic solutions are advanced directly on a Fourier grid.
- Dirichlet solutions are represented by an odd extension onto a periodic
  domain of length `2L`. Odd symmetry keeps the values at `x = 0` and `x = L`
  equal to zero.
- The optional `Cs` setting adds a subgrid-scale (SGS) viscosity model.

The implementation does not apply spectral dealiasing. Choose `N` and `dt`
carefully and check convergence when using the solver for quantitative work.

## Advanced use: stepping without storing every state

For long simulations, use the lower-level functions to process or save each
state as it is produced. This periodic example keeps only the current state:

```matlab
L = 64;
N = 128;
dt = 0.01;

x = L*(1:N)'/N;
state = sin(2*pi*x/L);
coeff = buildETDRK4('periodic',L,N,dt,1,1);

for step = 1:1000
    [~,flow] = stepKSPeriodic(coeff,L,N,state);
    state = flow(:,1);
end
```

For Dirichlet boundaries, the state passed to `stepKSDirichlet` contains only
the `N` interior values. Its output includes both zero endpoints:

```matlab
L = 64;
N = 128;
dt = 0.01;

xInterior = (1:N)'*L/(N+1);
state = sin(pi*xInterior/L);
coeff = buildETDRK4('dirichlet',L,N,dt,1,1);

for step = 1:1000
    [~,flow] = stepKSDirichlet(coeff,L,N,state);
    state = flow(2:end-1,1);
end
```

The first column of `flow` is the new solution. Columns two through five are
its first through fourth spatial derivatives. When `Cs` is supplied, column
six contains the SGS viscosity. Both step functions accept an optional final
`Cs` argument. The time-step size
is stored in `coeff` when `buildETDRK4` is called. Rebuild the coefficients
before changing the time-step size.

## Repository layout

| Path | Purpose |
| --- | --- |
| `src/solveKS.m` | Main interface for configuring and running a simulation. |
| `src/buildETDRK4.m` | Builds the spectral grid and ETDRK4 coefficients. |
| `src/stepKSPeriodic.m` | Advances a periodic solution by one time step. |
| `src/stepKSDirichlet.m` | Advances a Dirichlet solution by one time step. |
| `src/plotKSSpacetime.m` | Plots a solution with the original grayscale lit-surface style. |
| `examples/compare_boundaries.m` | Compares the two boundary conditions. |
