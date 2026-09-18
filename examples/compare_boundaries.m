addpath('../src');

config.L = 64;
config.N = 128;
config.dt = 0.1;
config.steps = 100;
config.initial = @(x) 0.5*sin(pi*x/config.L).*(1+0.3*sin(2*pi*x/config.L));

figure;
for j = 1:2
    if j == 1
        config.boundary = 'periodic';
    else
        config.boundary = 'dirichlet';
    end
    [t,x,u] = solveKS(config);
    subplot(1,2,j);
    imagesc(t,x,u);
    axis xy;
    xlabel('Time');
    ylabel('x');
    title(config.boundary);
end
