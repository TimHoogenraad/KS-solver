repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot,'src'));

% Long enough to show the initial transient and developed chaotic dynamics.
config.L = 128;
config.N = 256;
config.dt = 0.25;
config.steps = 2000;
config.initial = @(x) 0.5*sin(2*pi*x/config.L).*(1+0.3*sin(4*pi*x/config.L));

cases = struct( ...
    'boundary',{'periodic','dirichlet','dirichlet'}, ...
    'advection',{0,0,0.5}, ...
    'title',{'Periodic','Homogeneous Dirichlet', ...
             'Homogeneous Dirichlet, c = 0.5'});
solutions = cell(size(cases));
grids = cell(size(cases));

for j = 1:numel(cases)
    config.boundary = cases(j).boundary;
    config.advection = cases(j).advection;
    [t,grids{j},solutions{j}] = solveKS(config);
end

limit = max(cellfun(@(u) max(abs(u),[],'all'),solutions));
figure('Color','w','Units','inches','Position',[0 0 10 9.5]);
layout = tiledlayout(numel(cases),1,'TileSpacing','compact','Padding','compact');

for j = 1:numel(cases)
    ax = nexttile(layout);
    plotKSSpacetime(t,grids{j},solutions{j}, ...
        'Parent',ax, ...
        'DomainLength',config.L, ...
        'ColorLimit',limit, ...
        'Title',cases(j).title);
end

sgtitle(layout,'Kuramoto-Sivashinsky boundary and advection comparison');
