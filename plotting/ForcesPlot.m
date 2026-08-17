function ForcesPlot(Bend, Buoyancy, Shear, Out)

F_Visc = Bend.F_Visc;
F_Fail = Bend.F_Fail;
F_Bend = Bend.F_Bend;
F_Buoy = Buoyancy.F_Bouy;

ArcLength = Out.Slab.ArcLength;
[ArcLengthSubsumed, k] = SlabSubsumption(Out);

%%%-------------------------------------------------------------------------------%%%
%%% Plot the forces.
%%%-------------------------------------------------------------------------------%%%
font = 'Palatino Linotype';
fig = figure('DefaultTextFontName',font,'DefaultAxesFontName',font,...
    'DefaultAxesFontSize', 16);
theme(fig, "light")
set(gcf,'Color','w')
tiledlayout(2,1)

%%% Replace any nan values with the value at the base of the plate interface.
F_Horizontal = Shear.F_Horizontal;
i_nan = ~isfinite(F_Horizontal);
F_Horizontal(i_nan) = max(F_Horizontal);

nexttile
plot(1e-3*ArcLength, 1e-9*F_Fail, 'b', LineWidth = 3)
hold on
plot(1e-3*ArcLength, 1e-9*F_Horizontal, 'c', LineWidth = 3)
plot(1e-3*ArcLength, 1e-9*F_Buoy, 'm', LineWidth = 3)
y1 = round(1e-9*min(min([F_Buoy, F_Fail, F_Horizontal])) - 1);
y2 = round(1e-9*max(max([F_Buoy, F_Fail, F_Horizontal])) + 1);
plot(1e-3*ArcLength(k)*[1 1], [y1 y2], 'k--', LineWidth = 2)
hold off

ax = gca;
ax.YLim = [y1 y2];
xlabel('Arc Length (km)', 'Interpreter', 'latex')
ylabel('Force (GN/m)', 'Interpreter', 'latex')
legend('Bending', 'Shear Resistance',  'Bouyancy', 'Slab Subsumed', 'Interpreter',...
    'latex', Location = 'northwest')

nexttile
plot(1e-3*ArcLength, 1e-9*F_Fail, 'b', LineWidth = 3)
hold on
plot(1e-3*ArcLength, 1e-9*F_Visc, 'r', LineWidth = 3)
plot(1e-3*ArcLength, 1e-9*F_Bend, 'g', LineWidth = 3)

ax = gca;
ax.YScale = 'log';
xlabel('Arc Length (km)', 'Interpreter', 'latex')
ylabel('Bending Force (GN/m)', 'Interpreter', 'latex')

legend('Failure Envelope', 'Viscous', 'Viscous, H \& P (2019)', 'Interpreter',...
    'latex', Location = 'southeast')

axesH = findall(fig, "Type", "axes");
set(axesH, "TickLabelInterpreter", 'latex')