# KS solver

MATLAB ETDRK4 solver for the Kuramoto-Sivashinsky equation with a selectable
periodic or homogeneous Dirichlet boundary condition. The step functions are
copied from the latest `FullRunPeriodic_V2.m` and `FullRunDirichlet_V2.m`
workflows in `Documents/Research/KS_work/KS_cont`; their coefficient setup is
extracted into `buildETDRK4.m`. This repository contains solver code only.

## Run

Add `src` to the MATLAB path, then:

```matlab
config.boundary = 'dirichlet'; % or 'periodic'
config.L = 64;
config.N = 128;
config.dt = 0.1;
config.steps = 100;
config.initial = @(x) 0.5*sin(pi*x/config.L).*(1+0.3*sin(2*pi*x/config.L));
[t,x,u] = solveKS(config);
```

`u(:,j)` is the solution at `t(j)`. For periodic runs, `N` is the number of
grid points and `x` has `N` entries, including `L` but not `0`, matching the
source solver. For Dirichlet runs, `N` is the number of interior points and
`x` and `u` include the two zero-valued endpoints (`N+2` entries). A numeric
initial condition must have the same length as `x`.

Optional `config.Cs` enables the original SGS term. Optional `config.v2` and
`config.v4` set the linear coefficients, both `1` by default, as in the latest
long-run scripts. The nonlinear coefficient also follows `v2`, matching those
scripts. `config.steps` controls the number of time steps; all states are held
in memory, so keep it moderate or use the step functions directly for long runs.

Run `examples/compare_boundaries.m` from the `examples` directory for a plot.

## Scope

`solveKSV2.m` and `solveKSDirV2.m` preserve the long-run numerical methods,
including their optional SGS calculations. The original long-run scripts also
compute statistics, pQoI, and saved data; those are outside this solver repo.
The original scripts construct a dealias mask but their step functions do not
apply it. This package keeps that behavior rather than claiming dealiasing.
