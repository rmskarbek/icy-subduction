function Bend = BendingForce(Out)

%%% This function computes the slab bending force using three different methods.

%%%------------------------------------------------------------------------------------%%%
%%% OUTPUT. All of these variables are stored in a structure called 'Bend'.
%%%------------------------------------------------------------------------------------%%%
% dedt_ss       - Kinematically determined strain rate in the s-direction. 

% Sigma_ss_Visc - Normal stress in the s-direction, assuming that the slab is a viscous
%                 fluid with variable viscosity.

% Sigma_ss_Fail - Normal stress in the s-direction, when the viscous stresses are limited
%                 by the failure envelope.

% M_Visc        - Bending moment, assuming that the slab is a viscous fluid with variable
%                 viscosity.

% M_Fail        - Bending momemt when the viscous stresses are limited by the failure
%                 envelope.

% F_Visc        - Bending force, assuming that the slab is a viscous fluid with variable
%                 viscosity.

% F_Fail        - Bending force  when the viscous stresses are limited by the failure
%                 envelope.

% F_bend        - Bending force, assuming that the slab is a viscous fluid with variable
%                 viscosity, but using the incorrect calculation method from Howell & 
%                 Pappalardo (2019).

% S_C           - Failure envelopes for compression.

% S_T           - Failure envelopes for tension.

%%% NOTE: The bending forces (F_Visc, F_Fail, F_bend) are all horizontal forces resolved
%%% at the trench. That is, they are forces in the x-direction, at s = 0.


%%%------------------------------------------------------------------------------------%%%
%%% Geometry
%%%------------------------------------------------------------------------------------%%%
p = Out.p;
H = p.Geometry.SlabThick;               % thickness of conductive slab [m]
R_min = p.Geometry.CurveRadius;         % minimum radius of plate curvature [m]
v_plate = p.Geometry.PlateRate;         % plate convergence rate [m/s]

N = size(Out.Slab.Depth,1);

%%%------------------------------------------------------------------------------------%%%
%%% Simulation Output.
%%%------------------------------------------------------------------------------------%%%
Time = Out.TimeYrs;
Temp = Out.Temperature;
Eta = Out.Viscosity;

ArcLength = Out.Slab.ArcLength;
Depth = Out.Slab.Depth;
% Depth_init = Depth(:,1);
K_Center = Out.Slab.Curvature;
dKds = Out.Slab.dKds;

%%%------------------------------------------------------------------------------------%%%
%%% Bending force for a thin viscous plate.
%%%------------------------------------------------------------------------------------%%%
% Eta(:) = 1e21;

%%% Kinematic strain rate [1/s]
z = linspace(-H/2, H/2, N)';                    % depth in slab relative to centerline [m]
dedt_ss = -v_plate*repmat(z, 1, numel(Time)).*repmat(dKds', N, 1);               % [1/s]

%%% Compute the viscous fiber stress.
Sigma_ss_Visc = 4*dedt_ss.*Eta;                                                % [Pa]

%%% Compute the bending moment as a function of centerline distance using eq. (2) in 
%%% Buffett (2006).
M_Visc = nan(numel(Time),1);
for i = 1:numel(Time)
    M_Visc(i,1) = trapz(z, z.*Sigma_ss_Visc(:,i));                 % [Pa m^2 = N]
end

%%% Compute the bending force following eq. (20) in Buffett (2006).
F_Visc = cumtrapz(ArcLength, K_Center.*gradient(M_Visc)./gradient(ArcLength));  % [N/m]


%%%------------------------------------------------------------------------------------%%%
%%% Bending force using the failure envelope.
%%%------------------------------------------------------------------------------------%%%
mu_flag = 'temperature';
Sigma_ss_Fail = nan(N, numel(Time));
S_C = nan(N, numel(Time));
S_T = nan(N, numel(Time));

%%% Compute the failure envelope for each slab column, and define the fiber stress as the
%%% minimum of the viscous stress and failure stresses.
for i = 1:numel(Time)
    [S_Compression, ~, S_Tension, ~, ~] =...
        IceFailureEnvelope(Depth(:,i), Temp(:,i), mu_flag, dedt_ss(:,i));
    S_C(:,i) = S_Compression;
    S_T(:,i) = S_Tension;
    
    i_C = dedt_ss(:,i) > 0;
    i_T = ~i_C;
    Sigma_ss_Fail(i_C,i) = min([Sigma_ss_Visc(i_C,i), S_Compression(i_C)], [], 2);
    Sigma_ss_Fail(i_T,i) = max([Sigma_ss_Visc(i_T,i), S_Tension(i_T)], [], 2);
end

%%% Compute the bending moment as a function of centerline distance using eq. (2) in 
%%% Buffett (2006).
M_Fail = nan(numel(Time),1);
for i = 1:numel(Time)
    M_Fail(i,1) = trapz(z, z.*Sigma_ss_Fail(:,i));                 % [Pa m^2 = N]
end

%%% Compute the bending force following eq. (20) in Buffett (2006).
F_Fail = cumtrapz(ArcLength, K_Center.*gradient(M_Fail)./gradient(ArcLength));  % [N/m]


%%%------------------------------------------------------------------------------------%%%
%%% Bending force calculation from SubductionMain.n, Howell & Papallardo, 2019.
%%%------------------------------------------------------------------------------------%%%
F_Bend = zeros(numel(Time), 1);
for i = 1:numel(Time)
% Calculation of bending force requires reversal of eta vector
    eta_bend  = sort(Eta(:,i));
    F_bd      = zeros(N, 1);           % F_bend(depth)

% Initial bending force for full slab at lowest viscosity
    F_bd(1)   = (1/12)*v_plate*eta_bend(1);
    
% Calculate incremental contributions from higher viscosities
    for j = 2:N
        n = j-1;
        F_bd(j) = (1/12)*v_plate*(eta_bend(j) - eta_bend(j-1))*((N-n)/N)^3;
    end
    
% Incremental bending force scaled by amount of bending.
    dL = mean(diff(ArcLength));
    F_Bend(i) = sum(F_bd) * (dL/(2*R_min));
end

% Bending. For linear eta(z) this is off by a factor of 2 at s = 2*R_min. Seems odd.
F_Bend = cumsum(F_Bend);

%%%------------------------------------------------------------------------------------%%%
%%% Output
%%%------------------------------------------------------------------------------------%%%
Bend.dedt_ss = dedt_ss;
Bend.Sigma_ss_Visc = Sigma_ss_Visc;
Bend.Sigma_ss_Fail = Sigma_ss_Fail;
Bend.M_Visc = M_Visc;
Bend.M_Fail = M_Fail;
Bend.F_Visc = F_Visc;
Bend.F_Fail = F_Fail;
Bend.F_Bend = F_Bend;
Bend.S_C = S_C;
Bend.S_T = S_T;

% %%%------------------------------------------------------------------------------------%%%
% %%% Plot some failure enevlopes.
% %%%------------------------------------------------------------------------------------%%%
% font = 'Palatino Linotype';
% fig = figure('DefaultTextFontName',font,'DefaultAxesFontName',font,...
%     'DefaultAxesFontSize', 16);
% set(gcf,'Color','w')
% tiledlayout(2,2)
% 
% % z = Depth_init - Depth_init(end)/2;
% 
% %%%------------------------------------------------------------------------------------%%%
% nexttile
% i = 1;
% plot(1e-6*Sigma_ss_Visc(:,i), Depth_init, 'g', LineWidth = 3)
% 
% ax = gca;
% ax.YDir = "reverse";
% xlabel('Stress (MPa)', 'Interpreter', 'latex')
% ylabel('In-Slab Depth (m)', 'Interpreter', 'latex')
% title('Arc Length = 0 km', 'Interpreter', 'latex')
% 
% %%%------------------------------------------------------------------------------------%%%
% nexttile
% i_C = dedt_ss(:,i) > 0;
% i_T = ~i_C;
% 
% plot(1e-6*Sigma_ss_Visc(:,i), Depth_init, 'g', LineWidth = 3)
% hold on
% plot(1e-6*S_C(i_C,i),  Depth_init(i_C), 'm', LineWidth = 3)
% plot(1e-6*S_T(i_T,i),  Depth_init(i_T), 'c', LineWidth = 3)
% plot(1e-6*Sigma_ss_Fail(:,i), Depth_init, 'k--', LineWidth = 3);
% hold off
% 
% ax = gca;
% ax.YDir = "reverse";
% ax.XLim = [round(1e-6*min(S_T(:,i)) - 1) round(1e-6*max(S_C(:,i)) + 1)];
% xlabel('Stress (MPa)', 'Interpreter', 'latex')
% ylabel('In-Slab Depth (m)', 'Interpreter', 'latex')
% title('Arc Length = 0 km', 'Interpreter', 'latex')
% 
% %%%------------------------------------------------------------------------------------%%%
% nexttile
% i = 240;
% i_C = dedt_ss(:,i) > 0;
% i_T = ~i_C;
% 
% plot(1e-6*Sigma_ss_Visc(:,i), Depth_init, 'g', LineWidth = 3)
% hold on
% plot(1e-6*S_C(i_C,i),  Depth_init(i_C), 'm', LineWidth = 3)
% plot(1e-6*S_T(i_T,i),  Depth_init(i_T), 'c', LineWidth = 3)
% plot(1e-6*Sigma_ss_Fail(:,i), Depth_init, 'k--', LineWidth = 3);
% hold off
% 
% ax = gca;
% ax.YDir = "reverse";
% ax.XLim = [round(min(1e-6*Sigma_ss_Visc(:,i)) - 1) 2];
% xlabel('Stress (MPa)', 'Interpreter', 'latex')
% ylabel('In-Slab Depth (m)', 'Interpreter', 'latex')
% title('Arc Length = 5.2 km', 'Interpreter', 'latex')
% 
% legend('Viscous Stress', 'Brittle Failure', 'Ductile Failure', 'Composite Rheology',...
%     'Interpreter', 'latex', 'location', 'southwest')
% 
% %%%------------------------------------------------------------------------------------%%%
% nexttile
% [~, i] = min(abs(ArcLength - 8e3));
% i_C = dedt_ss(:,i) > 0;
% i_T = ~i_C;
% 
% plot(1e-6*Sigma_ss_Visc(:,i), Depth_init, 'g', LineWidth = 3)
% hold on
% plot(1e-6*S_C(i_C,i),  Depth_init(i_C), 'm', LineWidth = 3)
% plot(1e-6*S_T(i_T,i),  Depth_init(i_T), 'c', LineWidth = 3)
% plot(1e-6*Sigma_ss_Fail(:,i), Depth_init, 'k--', LineWidth = 3);
% hold off
% 
% ax = gca;
% ax.YDir = "reverse";
% xlabel('Stress (MPa)', 'Interpreter', 'latex')
% ylabel('In-Slab Depth (m)', 'Interpreter', 'latex')
% title('Arc Length = 8 km', 'Interpreter', 'latex')
% 
% axesH = findall(fig, "Type", "axes");
% set(axesH, "TickLabelInterpreter", 'latex')
% 
% 
% %%%------------------------------------------------------------------------------------%%%
% %%% Plot the bending forces.
% %%%------------------------------------------------------------------------------------%%%
% font = 'Palatino Linotype';
% fig = figure('DefaultTextFontName',font,'DefaultAxesFontName',font,...
%     'DefaultAxesFontSize', 16);
% set(gcf,'Color','w')
% tiledlayout(2,1)
% 
% nexttile
% plot(1e-3*ArcLength, 1e-9*F_Fail, 'b', LineWidth = 3)
% xlabel('Arc Length (km)', 'Interpreter', 'latex')
% ylabel('Bending Force (GN/m)', 'Interpreter', 'latex')
% 
% nexttile
% plot(1e-3*ArcLength, 1e-9*F_Fail, 'b', LineWidth = 3)
% hold on
% plot(1e-3*ArcLength, 1e-9*F_Visc, 'r', LineWidth = 3)
% plot(1e-3*ArcLength, 1e-9*F_Bend, 'g', LineWidth = 3)
% 
% ax = gca;
% ax.YScale = 'log';
% xlabel('Arc Length (km)', 'Interpreter', 'latex')
% ylabel('Bending Force (GN/m)', 'Interpreter', 'latex')
% 
% legend('Composite Rheology', 'Viscous', 'Viscous, H \& P (2019)', 'Interpreter', 'latex')
% 
% axesH = findall(fig, "Type", "axes");
% set(axesH, "TickLabelInterpreter", 'latex')