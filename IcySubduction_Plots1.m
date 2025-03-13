function IcySubduction_Plots1(Out)

%%%------------------------------------------------------------------------------------%%%
%%% Simulation Output.
%%%------------------------------------------------------------------------------------%%%
Time = Out.TimeYrs;
Depth = Out.Slab.Depth;
Depth_init = Depth(:,1);
ArcLength = Out.Slab.ArcLength;
Temp = Out.Temperature;
Phi = Out.Porosity;
Eta = Out.Viscosity;

%%%------------------------------------------------------------------------------------%%%
%%% Plots after Johnson2017, Figure 3.
%%%------------------------------------------------------------------------------------%%%
%%% First create a dummy figure to get the contour line for Temp = 259 K. This boundary
%%% can be used to infer where the slab is fully subsumed into the convecting layer.
% T_subsumed = Out.p.Temperature.BasalTemp - 0.01;
T_subsumed = 259;
f = figure;
M = contourf(repmat(1e-3*ArcLength', numel(Depth_init), 1), repmat(Depth_init, 1,...
    numel(Time)), Temp, T_subsumed*[1 1]);
X_subsumed = M(1,2:end)';
Y_subsumed = M(2,2:end)';
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
% subplot(3,1,1)
nexttile
contourf(repmat(1e-3*ArcLength', numel(Depth_init), 1), repmat(Depth_init, 1, numel(Time)),...
    Temp, EdgeColor="none")
hold on
plot(X_subsumed, Y_subsumed, 'g--', LineWidth=2)
hold off

c = colorbar;
c.TickLabelInterpreter = 'latex';
ax = gca;
% cmap_T = crameri('-lajolla');
cmap_T = flip(colormap('hot'));
colormap(ax, cmap_T);

ax.TickLabelInterpreter = 'latex';
ax.YDir = 'reverse';
title('Temperature (K)', 'Interpreter', 'latex')
xlabel('Arc Length (km)', 'Interpreter', 'latex')
ylabel('In-Slab Depth (m)', 'Interpreter', 'latex')

%%%------------------------------------------------------------------------------------%%%
%%% Viscosity
%%%------------------------------------------------------------------------------------%%%
nexttile
contourf(repmat(1e-3*ArcLength', numel(Depth_init), 1), repmat(Depth_init, 1, numel(Time)),...
    log10(Eta), EdgeColor="none")
hold on
plot(X_subsumed, Y_subsumed, 'g--', LineWidth=2)
hold off

c = colorbar;
c.TickLabelInterpreter = 'latex';
ax = gca;
colormap(ax, 'copper');
% colormap(ax, crameri('glasgow'));
hold on

%%% Plot a line where the plate first starts to unbend.
R_min = Out.p.Geometry.CurveRadius;
i = find(Out.Slab.ArcLength > R_min, 1);
plot(1e-3*ArcLength(i)*[1 1], [0 Depth_init(end)], 'r--', LineWidth = 2)
% plot(1e-3*ArcLength(364)*[1 1], [0 Depth_init(end)], 'r--', LineWidth = 2)

ax.YDir = 'reverse';
ax.TickLabelInterpreter = 'latex';
% XX = ax.XLim;
title('log(Viscosity) (Pa s)', 'Interpreter', 'latex')
xlabel('Arc Length (km)', 'Interpreter', 'latex')
ylabel('In-Slab Depth (m)', 'Interpreter', 'latex')


%%%------------------------------------------------------------------------------------%%%
%%% Porosity
%%%------------------------------------------------------------------------------------%%%
nexttile
contourf(repmat(1e-3*ArcLength', numel(Depth_init), 1), repmat(Depth_init, 1, numel(Time)),...
    Phi, EdgeColor="none")
hold on
plot(X_subsumed, Y_subsumed, 'g--', LineWidth=2)
hold off

c = colorbar;
c.TickLabelInterpreter = 'latex';
ax = gca;
cmap_P = colormap(ax, 'winter');
% cmap_P = crameri('roma');
colormap(ax, cmap_P);

ax.YDir = 'reverse';
ax.TickLabelInterpreter = 'latex';
% XX = ax.XLim;
title('Porosity', 'Interpreter', 'latex')
xlabel('Arc Length (km)', 'Interpreter', 'latex')
ylabel('In-Slab Depth (m)', 'Interpreter', 'latex')

% figure
% subplot(3,1,3)
% nexttile
% % plot(1e-3*real(Z_Bottom), 1e-3*imag(Z_Bottom), 'k', LineWidth=2)
% plot(1e-3*ArcLength, 1e-3*imag(Z_Bottom), 'k', LineWidth=2)
% hold on
% % plot(1e-3*real(Z_Top), 1e-3*imag(Z_Top), 'k', LineWidth=2)
% plot(1e-3*ArcLength, 1e-3*imag(Z_Top), 'k', LineWidth=2)
% hold off
% ax = gca;
% ax.TickLabelInterpreter = 'latex';
% ax.YDir = 'reverse';
% ax.XLim = XX;
% xlabel('Arc Length (km)', 'Interpreter', 'latex')
% ylabel('Shell Depth $y$-axis (km)', 'Interpreter', 'latex')
