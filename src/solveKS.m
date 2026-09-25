function [t,x,u] = solveKS(config)
% Run the latest ETDRK4 KS stepper with periodic or Dirichlet boundaries.
% Required fields: boundary, L, N, dt, steps, initial.
% Optional fields: Cs, v2, v4 (defaults: no SGS, 1, 1).
required = {'boundary','L','N','dt','steps','initial'};
for j = 1:numel(required)
    if ~isfield(config, required{j})
        error('Missing config.%s', required{j});
    end
end

boundary = lower(char(config.boundary));
if ~ismember(boundary, {'periodic','dirichlet'})
    error('config.boundary must be periodic or dirichlet');
end
if ~isscalar(config.L) || config.L <= 0 || ~isscalar(config.dt) || config.dt <= 0
    error('config.L and config.dt must be positive scalars');
end
if ~isscalar(config.steps) || config.steps < 0 || config.steps ~= floor(config.steps)
    error('config.steps must be a nonnegative integer');
end
if ~isscalar(config.N) || config.N < 2 || config.N ~= floor(config.N)
    error('config.N must be an integer of at least two');
end
if strcmp(boundary,'periodic') && mod(config.N,2) ~= 0
    error('Periodic config.N must be even');
end

v2 = 1;
v4 = 1;
if isfield(config,'v2'), v2 = config.v2; end
if isfield(config,'v4'), v4 = config.v4; end
useSGS = isfield(config,'Cs') && ~isempty(config.Cs);
if useSGS && (~isscalar(config.Cs) || config.Cs < 0)
    error('config.Cs must be a nonnegative scalar');
end

N = config.N;
if strcmp(boundary,'periodic')
    x = config.L*(1:N)'/N;
else
    x = (0:N+1)'*config.L/(N+1);
end
if isa(config.initial,'function_handle')
    initial = config.initial(x);
else
    initial = config.initial(:);
end
if numel(initial) ~= numel(x)
    error('config.initial must have one value per returned grid point');
end
initial = initial(:);
if strcmp(boundary,'dirichlet')
    initial([1,end]) = 0;
end

KS = buildETDRK4(boundary, config.L, N, config.dt, v2, v4);
t = (0:config.steps)*config.dt;
u = zeros(numel(x), config.steps+1);
u(:,1) = initial;
for j = 1:config.steps
    if strcmp(boundary,'periodic')
        if useSGS
            [~,flow] = stepKSPeriodic(KS,config.L,N,u(:,j),config.Cs);
        else
            [~,flow] = stepKSPeriodic(KS,config.L,N,u(:,j));
        end
    else
        if useSGS
            [~,flow] = stepKSDirichlet(KS,config.L,N,u(2:end-1,j),config.Cs);
        else
            [~,flow] = stepKSDirichlet(KS,config.L,N,u(2:end-1,j));
        end
    end
    u(:,j+1) = flow(:,1);
end
end
