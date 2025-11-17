function Johnson2017_Plots1(Out, Rho_Anom)

%%%------------------------------------------------------------------------------------%%%
%%% Simulation Output.
%%%------------------------------------------------------------------------------------%%%
Time = Out.TimeYrs;
Depth = Out.Slab.Depth;
Depth_init = 1e-3*Depth(:,1);
Temp = Out.Temperature;
Phi = Out.Porosity;

%%%------------------------------------------------------------------------------------%%%
%%% Plots after Johnson2017, Figure 3.
%%%------------------------------------------------------------------------------------%%%
%%% First create a dummy figure to get the contour line for Temp = 259 K. This boundary
%%% can be used to infer where the slab is fully subsumed into the convecting layer.
% T_subsumed = Out.p.Temperature.BasalTemp - 0.01;
T_subsumed = 259;
f = figure;
M = contourf(repmat(1e-3*Time', numel(Depth_init), 1), repmat(Depth_init, 1, numel(Time)),...
    Temp, T_subsumed*[1 1]);
X_subsumed = M(1,2:end)';
Y_subsumed = M(2,2:end)';
i_Tsub = X_subsumed == T_subsumed;
X_subsumed(i_Tsub) = nan;
Y_subsumed(i_Tsub) = nan;
close(f);

%%% Set up the main, three panelled figure.
font = 'Palatino Linotype';
figure('DefaultTextFontName',font,'DefaultAxesFontName',font,...
    'DefaultAxesFontSize', 16);
% tiledlayout(3, 1, 'padding', 'compact', 'TileSpacing', 'compact');
tiledlayout(3,1)

%%%------------------------------------------------------------------------------------%%%
%%% Temperature
%%%------------------------------------------------------------------------------------%%%
nexttile
% v = (100:20:260);
% contourf(repmat(1e-3*Time', numel(Depth_init), 1), repmat(Depth_init, 1, numel(Time)),...
%     Temp, v, EdgeColor="none")
% hold on
% plot(X_subsumed, Y_subsumed, 'g--', LineWidth=2)
% hold off

XX = [0 1e-3*Time(end)];
YY = [0 Depth_init(end)];
imagesc(XX, YY, Temp)
hold on
plot(X_subsumed, Y_subsumed, 'r--', 'LineWidth', 2)
hold off

c = colorbar;
c.TickLabelInterpreter = 'latex';
ax = gca;
% cmap_T = crameri('-lajolla');
% cmap_T = flip(colormap('hot'));
% colormap(ax, cmap_T);

ax.TickLabelInterpreter = 'latex';
ax.XLim = XX;
% ax.YDir = 'reverse';
ax.YTick = (0:1.5:6);
title('Temperature (K)', 'Interpreter', 'latex')
xlabel('Time (kyr)', 'Interpreter', 'latex')
ylabel('In-Slab Depth (m)', 'Interpreter', 'latex')


%%%------------------------------------------------------------------------------------%%%
%%% Porosity
%%%------------------------------------------------------------------------------------%%%
nexttile
% contourf(repmat(1e-3*Time', numel(Depth_init), 1), repmat(Depth_init, 1, numel(Time)),...
%     Phi, EdgeColor="none")
% hold on
% plot(X_subsumed, Y_subsumed, 'g--', LineWidth=2)
% hold off

XX = [0 1e-3*Time(end)];
YY = [0 Depth_init(end)];
imagesc(XX, YY, Phi)
hold on
plot(X_subsumed, Y_subsumed, 'r--', 'LineWidth', 2)
hold off

c = colorbar;
c.TickLabelInterpreter = 'latex';
ax = gca;
% cmap_P = colormap(ax, 'winter');
% cmap_P = crameri('roma');
colormap(ax, jet);

% ax.YDir = 'reverse';
ax.XLim = XX;
ax.TickLabelInterpreter = 'latex';
ax.YTick = (0:1.5:6);
title('Porosity', 'Interpreter', 'latex')
xlabel('Time (kyr)', 'Interpreter', 'latex')
ylabel('In-Slab Depth (m)', 'Interpreter', 'latex')

%%%------------------------------------------------------------------------------------%%%
%%% Geometry
%%%------------------------------------------------------------------------------------%%%
Z_Bottom = Out.Slab.Bottom;
Z_Top = Out.Slab.Top;

nexttile
plot(1e-3*Time, 1e-3*imag(Z_Bottom), 'r', LineWidth=2)
hold on
plot(1e-3*Time, 1e-3*imag(Z_Top), 'b', LineWidth=2)
hold off

ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.YDir = 'reverse';
ax.YTick = (0:5:20);
ax.XLim = [0 1e-3*Time(end)];
xlabel('Time (kyr)', 'Interpreter', 'latex')
ylabel('Shell Depth $y$-axis (km)', 'Interpreter', 'latex')

%%%------------------------------------------------------------------------------------%%%
%%% Density anomaly
%%%------------------------------------------------------------------------------------%%%
figure('DefaultTextFontName',font,'DefaultAxesFontName',font, 'DefaultAxesFontSize', 16);

plot(1e-3*Time, Rho_Anom, 'k', 'LineWidth', 2);
hold on
plot(1e-3*[0 Time(end)], [0 0], 'k')
hold off

ax = gca;
ax.TickLabelInterpreter = 'latex';
xlabel('Time (kyr)', 'Interpreter', 'latex')
ylabel('Average Density Anomaly (kg/m$^3$)', 'Interpreter', 'latex')