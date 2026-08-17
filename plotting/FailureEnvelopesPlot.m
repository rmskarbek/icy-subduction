function FailureEnvelopesPlot(Bend, Out)

dedt_ss = Bend.dedt_ss;
Sigma_ss_Visc = Bend.Sigma_ss_Visc;
Sigma_ss_Fail = Bend.Sigma_ss_Fail;

S_C = Bend.S_C;
S_T = Bend.S_T;

Depth = Out.Slab.Depth;
Depth_init = Depth(:,1);
Depth_init = Depth_init - Depth_init(end)/2;

%%%-------------------------------------------------------------------------------%%%
%%% Initialize the figure.
%%%-------------------------------------------------------------------------------%%%
font = 'Palatino Linotype';
fig = figure('DefaultTextFontName',font,'DefaultAxesFontName',font,...
    'DefaultAxesFontSize', 16);
theme(fig, "light")
set(gcf,'Color','w')
tiledlayout(2,2)

%%%-------------------------------------------------------------------------------%%%
%%% Plot just the viscous stress profile at the trench to show how large it is.
%%%-------------------------------------------------------------------------------%%%
nexttile
i = 1;
plot(1e-6*Sigma_ss_Visc(:,i), Depth_init, 'g', LineWidth = 3)

ax = gca;
ax.YDir = "reverse";
xlabel('Stress (MPa)', 'Interpreter', 'latex')
ylabel('In-Slab Depth (m)', 'Interpreter', 'latex')
title('Arc Length = 0 km', 'Interpreter', 'latex')

%%%-------------------------------------------------------------------------------%%%
%%% Plot the failure envelope stresses at the trench to show how smaller they are 
%%% compared to the viscous stress.
%%%-------------------------------------------------------------------------------%%%
nexttile
i_C = dedt_ss(:,i) > 0;
i_T = ~i_C;

plot(1e-6*Sigma_ss_Visc(:,i), Depth_init, 'g', LineWidth = 3)
hold on
plot(1e-6*S_C(i_C,i),  Depth_init(i_C), 'm', LineWidth = 3)
plot(1e-6*S_T(i_T,i),  Depth_init(i_T), 'c', LineWidth = 3)
plot(1e-6*Sigma_ss_Fail(:,i), Depth_init, 'k--', LineWidth = 3);
hold off

ax = gca;
ax.YDir = "reverse";
ax.XLim = [round(1e-6*min(S_T(i_T,i)) - 1) round(1e-6*max(S_C(i_C,i)) + 1)];
xlabel('Stress (MPa)', 'Interpreter', 'latex')
ylabel('In-Slab Depth (m)', 'Interpreter', 'latex')
title('Arc Length = 0 km', 'Interpreter', 'latex')

%%%-------------------------------------------------------------------------------%%%
%%% Plot the failure envelop stresses at the arc length where the plate starts to 
%%% unbend.
%%%-------------------------------------------------------------------------------%%%
R_min = Out.p.Geometry.CurveRadius;
i = find(Out.Slab.ArcLength > R_min, 1);
i_C = dedt_ss(:,i) > 0;
i_T = ~i_C;

nexttile
plot(1e-6*Sigma_ss_Visc(:,i), Depth_init, 'g', LineWidth = 3)
hold on
plot(1e-6*S_C(i_C,i),  Depth_init(i_C), 'm', LineWidth = 3)
plot(1e-6*S_T(i_T,i),  Depth_init(i_T), 'c', LineWidth = 3)
plot(1e-6*Sigma_ss_Fail(:,i), Depth_init, 'k--', LineWidth = 3);
hold off

ax = gca;
ax.YDir = "reverse";
ax.XLim = [round(1e-6*min(S_T(i_T,i)) - 1) round(1e-6*max(S_C(i_C,i)) + 1)];
xlabel('Stress (MPa)', 'Interpreter', 'latex')
ylabel('In-Slab Depth (m)', 'Interpreter', 'latex')
title(sprintf('Arc Length = %0.3f km', 1e-3*Out.Slab.ArcLength(i)), 'Interpreter',...
    'latex')

legend('Viscous Stress', 'Brittle Failure', 'Ductile Failure', 'Composite Rheology',...
    'Interpreter', 'latex', 'location', 'southwest')

%%%-------------------------------------------------------------------------------%%%
%%% Plot the failure envelope stresses where the plate is completely subsumed, or at 
%%% the maximum arc length of the slab.
%%%-------------------------------------------------------------------------------%%%
[~, ~, ~, k, ~] = SlabSubsumption(Out);
nexttile
i_C = dedt_ss(:,k) > 0;
i_T = ~i_C;

plot(1e-6*Sigma_ss_Visc(:,k), Depth_init, 'g', LineWidth = 3)
hold on
plot(1e-6*S_C(i_C,k),  Depth_init(i_C), 'm', LineWidth = 3)
plot(1e-6*S_T(i_T,k),  Depth_init(i_T), 'c', LineWidth = 3)
plot(1e-6*Sigma_ss_Fail(:,k), Depth_init, 'k--', LineWidth = 3);
hold off

ax = gca;
ax.YDir = "reverse";
% ax.XLim = [round(1e-6*min(S_T(:,k)) - 1) round(1e-6*max(S_C(:,k)) + 1)];
xlabel('Stress (MPa)', 'Interpreter', 'latex')
ylabel('In-Slab Depth (m)', 'Interpreter', 'latex')
title(sprintf('Arc Length = %0.3f km', 1e-3*Out.Slab.ArcLength(k)), 'Interpreter', 'latex')

axesH = findall(fig, "Type", "axes");
set(axesH, "TickLabelInterpreter", 'latex')