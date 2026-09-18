# Kuramoto-Sivashinsky solver

This MATLAB repository solves the one-dimensional Kuramoto-Sivashinsky (KS)
equation with either periodic or homogeneous Dirichlet boundary conditions.
Set `config.boundary` to choose between them. The code uses the Fourier
spectral ETDRK4 methods from the newest periodic and Dirichlet long-run
workflows in `KS_cont`.

The package contains the solver only. It does not include the original
long-run analysis, pQoI or refinement code, figures, or generated `.mat` data.

## Repository contents

| File | Purpose |
| --- | --- |
| `src/solveKS.m` | Main entry point. Checks settings, builds the grid, advances the solution, and returns the full history. |
| `src/buildETDRK4.m` | Precomputes Fourier wave numbers and ETDRK4 coefficients for the selected boundary condition. |
| `src/stepKSPeriodic.m` | One periodic ETDRK4 step, with an optional subgrid-scale (SGS) term. |
| `src/stepKSDirichlet.m` | One Dirichlet ETDRK4 step, using an odd extension on a doubled periodic domain; optional SGS term. |
| `examples/compare_boundaries.m` | Runs and plots a short case with each boundary condition. |
| `.gitignore` | Excludes generated MATLAB data, figures, and autosave files. |

`solveKS` is the normal way to run a case. The `stepKS*` functions are useful
when writing a long-running driver that manages its own output and avoids
keeping every state in memory.

## Requirements

MATLAB with `fft` and `ifft`. No additional toolbox is used by the solver.
From MATLAB, add the `src` directory to the path:

```matlab
addpath('/path/to/ks-solver/src')
```

## Quick start

```matlab
config.boundary = 'periodic'; % change to 'dirichlet' for zero endpoints
config.L = 64;
config.N = 128;
config.dt = 0.01;
config.steps = 100;
config.initial = @(x) 0.5*sin(2*pi*x/config.L) ...
    .* (1 + 0.3*sin(4*pi*x/config.L));

[t,x,u] = solveKS(config);
imagesc(t,x,u); axis xy; xlabel('Time'); ylabel('x');
```

To run the included comparison from the repository root:

```matlab
run('examples/compare_boundaries.m')
```

The example finds `src` relative to its own location, so it also works when
launched from another MATLAB working directory with an absolute script path.

## Settings

| Field | Meaning |
| --- | --- |
| `boundary` | Required: `'periodic'` or `'dirichlet'`. |
| `L` | Required: domain length, a positive number. |
| `N` | Required: periodic grid points (must be even), or Dirichlet **interior** grid points. |
| `dt` | Required: positive time-step size. |
| `steps` | Required: nonnegative integer number of time steps. |
| `initial` | Required: function of `x`, or a vector with one value per returned grid point. |
| `Cs` | Optional: nonnegative SGS coefficient. Omit it for a DNS-style run. |
| `v2`, `v4` | Optional: linear coefficients; both default to `1`, matching the latest long-run scripts. |

The returned `t` is a row vector of length `steps+1`; `u(:,j)` is the solution
at `t(j)`. The returned `x` is a column vector. With periodic boundaries, it
has `N` points at `L/N, 2L/N, ..., L`; it includes `L` but not `0`, following
the source implementation. With Dirichlet boundaries, it has `N+2` points
including `0` and `L`. The two endpoint values are set to zero, including at
the initial time. A numeric `initial` vector must match this returned grid.

`solveKS` stores all states in RAM. The memory for `u` grows with both `N`
and `steps`. For a long simulation, use the lower-level step functions and
save or analyze states incrementally.

For example, this advances a periodic state without building a history array:

```matlab
L = 64; N = 128; dt = 0.01;
x = L*(1:N)'/N;
state = sin(2*pi*x/L);
coeff = buildETDRK4('periodic',L,N,dt,1,1);
for j = 1:1000
    [~,flow] = stepKSPeriodic(coeff,L,N,dt,state);
    state = flow(:,1);
end
```

For Dirichlet, call `buildETDRK4('dirichlet',L,N,dt,1,1)` and
`stepKSDirichlet(coeff,L,N,dt,state)`, where `state` has `N` interior values.
For the next step, set `state = flow(2:end-1,1)` to drop the returned endpoints.
Both step functions accept an optional final `Cs` argument. The `h` argument
in a step call must equal the `dt` used to build `coeff`; the original step
functions take it for interface compatibility, while the precomputed
coefficients determine the actual step size.

## Numerical method and provenance

Both branches use ETDRK4 with a Fourier spectral spatial representation.
Periodic runs advance the field directly on a periodic grid. Dirichlet runs
oddly extend the interior field onto a domain of length `2L`, advance that
field, and enforce zero-valued endpoints after each step. The optional SGS
term is retained from the source functions.

`stepKSPeriodic.m` and `stepKSDirichlet.m` are the functions formerly named
`solveKSV2.m` and `solveKSDirV2.m` in `KS_cont/functions`. Their numerical
bodies were preserved; only the function and file names changed. The ETDRK4
coefficient formulas came from `FullRunPeriodic_V2.m` and
`FullRunDirichlet_V2.m`. Those source scripts create a dealiasing mask, but
their step functions never apply it. This package therefore does not claim
dealiasing.

The lower-level step functions return `[x,flow]`. The first column of `flow`
is the new field; the next four are spatial derivatives of orders one through
four. When `Cs` is supplied, column six is the SGS viscosity. For periodic
runs, `flow` has `N` rows. For Dirichlet runs, it has `N+2` rows with endpoints.
The high-level `solveKS` function returns only the field history.

Short periodic and Dirichlet runs, both with and without SGS, were checked
with MATLAB R2026a. The checks covered output dimensions, finite values,
and zero Dirichlet endpoints.
