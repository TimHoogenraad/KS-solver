function ax = plotKSSpacetime(t,x,u,varargin)
% Plot a KS solution using the original grayscale lit-surface style.
% Name-value options:
%   Parent       Axes to plot into. A new figure is created by default.
%   DomainLength Physical domain length. Defaults to max(x).
%   ColorLimit   Symmetric color limit. Defaults to max(abs(u)).
%   Title        Axes title.

parser = inputParser;
addParameter(parser,'Parent',[],@(value) isempty(value) || isgraphics(value,'axes'));
addParameter(parser,'DomainLength',max(x),@(value) isscalar(value) && value > 0);
addParameter(parser,'ColorLimit',[],@(value) isempty(value) || (isscalar(value) && value > 0));
addParameter(parser,'Title','',@(value) ischar(value) || isstring(value));
parse(parser,varargin{:});

if size(u,1) ~= numel(x) || size(u,2) ~= numel(t)
    error('u must have numel(x) rows and numel(t) columns');
end

ax = parser.Results.Parent;
if isempty(ax)
    figure('Color','w','Units','inches','Position',[0 0 6.5 3.0]);
    ax = axes;
end

domainLength = parser.Results.DomainLength;
colorLimit = parser.Results.ColorLimit;
if isempty(colorLimit)
    colorLimit = max(abs(u),[],'all');
    if colorLimit == 0
        colorLimit = 1;
    end
end

surf(ax,t,x,u);
shading(ax,'interp');
lighting(ax,'gouraud');
axis(ax,'tight');
view(ax,2);
zlim(ax,[-50 500]);
clim(ax,[-colorLimit colorLimit]);
colormap(ax,gray(256));
material(ax,[0.30 0.60 0.60 40.00 1.00]);
light(ax,'Position',[0,0,2],'Style','infinite');

xlabel(ax,'Time');
ylabel(ax,'x/L');
title(ax,parser.Results.Title);
yticks(ax,linspace(0,domainLength,5));
yticklabels(ax,{'0','1/4','1/2','3/4','1'});
ylim(ax,[0 domainLength]);
set(ax,'Layer','top','TickDir','out','FontSize',11);
end
