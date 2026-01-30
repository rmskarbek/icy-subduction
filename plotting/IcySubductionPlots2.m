function IcySubductionPlots2(Out)

%%% Get the geometry of the slab.
p = Out.p;
H = p.Geometry.SlabThick;
H_shell = p.Geometry.ShellThick;

%%% Get the slab geometry mesh grid.
[x_Slab, y_Slab] = SlabGrid(Out);

%%% Convert units to kilometers.
x_Slab = 1e-3*x_Slab;
y_Slab = 1e-3*y_Slab;
H = 1e-3*H;
H_shell = 1e-3*H_shell;

%%%-------------------------------------------------------------------------------%%%
%%% Temperature
%%%-------------------------------------------------------------------------------%%%
%%% Set up the figure.
font = 'Palatino Linotype';
f = figure('DefaultTextFontName',font,'DefaultAxesFontName',font,...
    'DefaultAxesFontSize', 16);
% theme(f, "dark")
theme(f, "light")
% set(gcf,'Color','w')

%%% Set the color map.
cmap_T = crameri('lajolla');

%%% Plot the convective layer.
x1 = 0;
x2 = ceil(max(max(x_Slab)));

area([x1 x2], H_shell*[1 1], FaceColor = cmap_T(end,:));
hold on

%%% Plot the conductive layer.
% area([x1 x2], H*[1 1], FaceColor = cmap_T(1,:));
nz = size(x_Slab,1);
ns = size(x_Slab,2);
x_Cond = [x_Slab(1,1)*ones(nz,1), x_Slab(1,round(ns/2))*ones(nz,1),...
    x_Slab(1,end)*ones(nz,1)];
y_Cond = y_Slab(:,1)*[1 1 1];

TemperatureCond = Out.Temperature(:,1)*[1 1 1];
contourf(x_Cond, y_Cond, TemperatureCond, EdgeColor="none")
plot([1 x2], H*[1 1], 'c', LineWidth=1)

%%% Plot the slab temperature field. Round the temperature to the nearest integer to
%%% remove unwanted countour line noise.
contourf(x_Slab, y_Slab, round(Out.Temperature), EdgeColor="none")
c = colorbar;
c.TickLabelInterpreter = 'latex';
ax = gca;
colormap(ax, cmap_T);
ax.YDir = "reverse";
axis equal

%%% Plot outlines of the slab so that it is more visible.
plot(x_Slab(1,:), y_Slab(1,:), 'c', LineWidth=1)
plot(x_Slab(end,:), y_Slab(end,:), 'c', LineWidth=1)
plot(x_Slab(:,end), y_Slab(:,end), 'c', LineWidth=1)

%%% Plot subsumption location.
[~, i_Subsumed] = SlabSubsumption(Out);
x_Subsumed = 1e-3*real([Out.Slab.Bottom(i_Subsumed), Out.Slab.Top(i_Subsumed)]);
y_Subsumed = 1e-3*imag([Out.Slab.Bottom(i_Subsumed), Out.Slab.Top(i_Subsumed)]);
plot(x_Subsumed, y_Subsumed, 'b')
hold off

%%% Title and labels.
title('Temperature (K)', 'Interpreter', 'latex')
xlabel('Distance (km)', 'Interpreter', 'latex')
ylabel('Depth (km)', 'Interpreter', 'latex')


%%%-------------------------------------------------------------------------------%%%
%%% Porosity
%%%-------------------------------------------------------------------------------%%%
%%% Set up the figure.
font = 'Palatino Linotype';
f = figure('DefaultTextFontName',font,'DefaultAxesFontName',font,...
    'DefaultAxesFontSize', 16);
% theme(f, "dark")
theme(f, "light")
% set(gcf,'Color','w')

%%% Set the color map.
cmap_P = crameri('roma');

%%% Plot the convective layer.
x1 = 0;
x2 = max(max(x_Slab));

area([x1 x2], H_shell*[1 1], FaceColor = cmap_P(1,:));
hold on

%%% Plot the conductive layer.
% area([x1 x2], H*[1 1], FaceColor = cmap_P(end,:));
PorosityCond = Out.Porosity(:,1)*[1 1 1];
contourf(x_Cond, y_Cond, PorosityCond, EdgeColor="none")
plot([1 x2], H*[1 1], 'w', LineWidth=1)

%%% Plot the slab temperature field. Round the temperature to the nearest integer to
%%% remove unwanted countour line noise.
contourf(x_Slab, y_Slab, Out.Porosity, EdgeColor="none")
c = colorbar;
c.TickLabelInterpreter = 'latex';
ax = gca;
colormap(ax, cmap_P);
ax.YDir = "reverse";
axis equal

%%% Plot outlines of the slab so that it is more visible.
plot(x_Slab(1,:), y_Slab(1,:), 'w', LineWidth=1)
plot(x_Slab(end,:), y_Slab(end,:), 'w', LineWidth=1)
plot(x_Slab(:,end), y_Slab(:,end), 'w', LineWidth=1)
hold off

%%% Title and labels.
title('Porosity', 'Interpreter', 'latex')
xlabel('Distance (km)', 'Interpreter', 'latex')
ylabel('Depth (km)', 'Interpreter', 'latex')