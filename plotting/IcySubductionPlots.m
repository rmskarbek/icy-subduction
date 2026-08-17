function IcySubductionPlots(Out)

%%% This function will make contour plots of the (1) temperature, (2) porosity, and
%%% (3) bulk density fields from a subduction simulation.

%%%-------------------------------------------------------------------------------%%%
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
set(gcf,'Color','w')

%%% Set the color map.
cmap_T = crameri('lajolla');

%%% Plot the convective layer.
x1 = 0;
x2 = max(max(x_Slab));

area([x1 x2], H_shell*[1 1], FaceColor = cmap_T(end,:));
hold on

%%% Plot the conductive layer.
nz = size(x_Slab,1);
ns = size(x_Slab,2);
x_Cond = [x_Slab(1,1)*ones(nz,1), x_Slab(1,round(ns/2))*ones(nz,1),...
    x_Slab(1,end)*ones(nz,1)];
y_Cond = y_Slab(:,1)*[1 1 1];

[~, ~, ~, i_Subsumed, T_subsumed] = SlabSubsumption(Out);
temp_levels = [100:20:240, T_subsumed];
TemperatureCond = Out.Temperature(:,1)*[1 1 1];
contourf(x_Cond, y_Cond, TemperatureCond, temp_levels, EdgeColor="none")
plot([1 x2], H*[1 1], 'c', LineWidth=1)

%%% Plot the slab temperature field.
contourf(x_Slab, y_Slab, Out.Temperature, temp_levels, EdgeColor="none")
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
set(gcf,'Color','w')

%%% Set the color map.
cmap_P = crameri('roma');

%%% Plot the convective layer.
x1 = 0;
x2 = max(max(x_Slab));

area([x1 x2], H_shell*[1 1], FaceColor = cmap_P(1,:));
hold on

%%% Plot the conductive layer.
PorosityCond = Out.Porosity(:,1)*[1 1 1];
contourf(x_Cond, y_Cond, PorosityCond, EdgeColor="none")
plot([1 x2], H*[1 1], 'w', LineWidth=1)

%%% Plot the slab porosity field.
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


%%%-------------------------------------------------------------------------------%%%
%%% Bulk Density
%%%-------------------------------------------------------------------------------%%%
%%% Set up the figure.
font = 'Palatino Linotype';
f = figure('DefaultTextFontName',font,'DefaultAxesFontName',font,...
    'DefaultAxesFontSize', 16);
% theme(f, "dark")
theme(f, "light")
set(gcf,'Color','w')

%%% Set the color map.
cmap_D = crameri('-batlow');
% cmap_D = crameri('-acton');
% cmap_D = cmap_D(64:end,:);

%%% Determine the color for the convecting ice.
rho_max = max(max(Out.Density));
rho_min = min(min(Out.Density));
rho_conv = IceDensity(Out.p.Temperature.BasalTemp, 'density');
i_conv = ceil(((rho_conv - rho_min)/(rho_max - rho_min))*size(cmap_D,1));

%%% Plot the convective layer.
x1 = 0;
x2 = max(max(x_Slab));
area([x1 x2], H_shell*[1 1], FaceColor = cmap_D(i_conv,:));
hold on

%%% Define density contour levels to pick out the convective ice density.
% dens_levels = [rho_min, (840:10:910), rho_conv, 920, 930];
dens_levels = [rho_min, (840:10:910), 918, 920, 930];

%%% Plot the conductive layer.
DensityCond = Out.Density(:,1)*[1 1 1];
% contourf(x_Cond, y_Cond, DensityCond, EdgeColor="none")
contourf(x_Cond, y_Cond, DensityCond, dens_levels, EdgeColor="none")
plot([1 x2], H*[1 1], 'c', LineWidth=1)

%%% Plot the slab bulk density field.
% contourf(x_Slab, y_Slab, Out.Density, EdgeColor="none");
contourf(x_Slab, y_Slab, Out.Density, dens_levels, EdgeColor="none");
% contour(x_Slab, y_Slab, Out.Density, 920*[1 1], EdgeColor="y")
c = colorbar;
c.TickLabelInterpreter = 'latex';
ax = gca;
colormap(ax, cmap_D);
ax.YDir = "reverse";
axis equal

%%% Plot outlines of the slab so that it is more visible.
plot(x_Slab(1,:), y_Slab(1,:), 'c', LineWidth=1)
plot(x_Slab(end,:), y_Slab(end,:), 'c', LineWidth=1)
plot(x_Slab(:,end), y_Slab(:,end), 'c', LineWidth=1)
hold off

%%% Title and labels.
title('Bulk Density (kg/m$^3$)', 'Interpreter', 'latex')
xlabel('Distance (km)', 'Interpreter', 'latex')
ylabel('Depth (km)', 'Interpreter', 'latex')
