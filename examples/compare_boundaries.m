repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot,'src'));

% Long enough to show the initial transient and developed chaotic dynamics.
config.L = 128;
config.N = 256;
config.dt = 0.25;
config.steps = 2000;
config.initial = @(x) 0.5*sin(2*pi*x/config.L).*(1+0.3*sin(4*pi*x/config.L));

boundaries = {'periodic','dirichlet'};
titles = {'Periodic','Homogeneous Dirichlet'};
solutions = cell(1,2);
grids = cell(1,2);

for j = 1:2
    config.boundary = boundaries{j};
    [t,grids{j},solutions{j}] = solveKS(config);
end

limit = max(cellfun(@(u) max(abs(u),[],'all'),solutions));
figure('Color','w','Units','inches','Position',[0 0 10 6.5]);
layout = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

for j = 1:2
    ax = nexttile(layout);
    plotKSSpacetime(t,grids{j},solutions{j}, ...
        'Parent',ax, ...
        'DomainLength',config.L, ...
        'ColorLimit',limit, ...
        'Title',titles{j});
end

sgtitle(layout,'Kuramoto-Sivashinsky transient and chaotic dynamics');
