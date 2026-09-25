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
    surf(ax,t,grids{j}/config.L,solutions{j},'EdgeColor','none');
    shading(ax,'interp');
    view(ax,2);
    axis(ax,'tight');
    clim(ax,[-limit limit]);
    colormap(ax,gray(256));
    material(ax,[0.30 0.60 0.60 40.00 1.00]);
    light(ax,'Position',[0,0,2],'Style','infinite');
    title(ax,titles{j});
    ylabel(ax,'x/L');
    yticks(ax,[0 0.25 0.5 0.75 1]);
    set(ax,'Layer','top','TickDir','out','FontSize',11);
end

xlabel(nexttile(layout,2),'Time');
sgtitle(layout,'Kuramoto-Sivashinsky transient and chaotic dynamics');
