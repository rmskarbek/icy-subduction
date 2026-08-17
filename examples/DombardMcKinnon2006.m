function [Depth, S] = DombardMcKinnon2006
%%% This function will compute a failure envelope for Europa, following Dombard &
%%% McKinnon (2006) and reproducing figure 2 in that paper.

%%% Some constants from Dombard & McKinnon.
D = 5;                                                      % shell thickness [km]
g = 1.31;                                                   % Europa gravity [m / s^2]
R = 8.314;                                                  % gas constant [J / (mol K]
rho_ice = 950;                                              % ice density [kg/m^3]

%%% Depth [km].
Depth = linspace(0, D, 1000)';

%%% Temperature profile.
T_s = 80;                                                   % Surface temperature [K]
dTdz = 18;                                                  % thermal gradient [K/km]
TempK = T_s + dTdz*Depth;                                   % [K]

%%% Lithostatic stress [MPa]
sigma_L = 1e-3*g*rho_ice*Depth;

%%% Friction failure in terms of principal stresses, assuming that sigma_3 is equal to
%%% the lithostatic stress. This also assumes that failure surfaces are on optimally 
%%% oriented planes. I.E., planes at an angle phi = atan(mu) measured from sigma_3 in 
%%% the direction of sigma_1. So the orientation of the failures surfaces are 
%%% determined by the friction coefficient.
mu_0 = 0.55;                                           % eq (1a) in Dombard & McKinnon
phi = atan(mu_0);
C = 1;                                                 % Cohesion [MPa]
sigma_1 = (2*C + (cos(phi) + mu_0.*(sin(phi) + 1)).*sigma_L)...
    ./(cos(phi) + mu_0.*(sin(phi) - 1));
Tau_F = sigma_1 - sigma_L;

%%% Ice flow laws, paramters for Regimes A, B, and C are from from Dombard & McKinnon 
%%% (2001), Table 1. These regimes represent grainsize insensitive disclocation creep 
%%% mechanisms. See Durham & Stern (2001) for references. See text in D & M (2006) in 
%%% the paragraph after eq (1b). Symbols here are after Dombard & McKinnon (2001).
e = 1e-15;                                             % strain rate [1/s]

%%% Regime A, 240 - 258 K.
Q_A = 91e3;                                             % [J / mol]
n_A = 4;                                                % stress exponent
A_A = 10^11.8;                                          % [1 / (s MPa^n)]
A_A = A_A*(3^(1/2)/2)^(n_A + 1);
Tau_A = ((e/A_A)*exp(Q_A./(R*TempK))).^(1/n_A);         % [MPa]

%%% Regime B.
Q_B = 61e3;                                             % [J / mol]
n_B = 4;                                                % stress exponent
A_B = 10^5.1;                                           % [1 / (s MPa^n)]
A_B = A_B*(3^(1/2)/2)^(n_B + 1);
Tau_B = ((e/A_B)*exp(Q_B./(R*TempK))).^(1/n_B);         % [MPa]

%%% Regime C.
% Q_C = 36e3;                                             % [J / mol]
% n_C = 4.7;
% A_C = 10^(-2.8);                                        % [1 / (s MPa^n)]
Q_C = 39e3;                                             % [J / mol]
n_C = 6;                                                % stress exponent
A_C = 10^(-3.8);                                        % [1 / (s MPa^n)]
A_C = A_C*(3^(1/2)/2)^(n_C + 1);
Tau_C = ((e/A_C)*exp(Q_C./(R*TempK))).^(1/n_C);         % [MPa]

%%% Grainsize sensitive grain boundary sliding from Goldsby & Kohlstedt (2001). Here,
%%% using parameter values listed in Table in in Dombard & McKinnon (2001).
d = 1e-3;                                               % grainsize [m]
Q_G = 49e3;                                             % [J / mol]
m_G = 1.4;                                              % grainsize exponent
n_G = 1.8;                                              % stress exponent
A_G = 3.9e-3;                                           % [1 / (s m^m MPa^n)]
A_G = A_G*(3^(1/2)/2)^(n_G + 1);
Tau_G = ((e*d^m_G/A_G)*exp(Q_G./(R*TempK))).^(1/n_G);         % [MPa]

%%% Find the failure envelope.
S = [Tau_F, Tau_A, Tau_B, Tau_C, Tau_G];
S = min(S,[],2);

%%% Plot the failure envelope. Compare with dashed line for 'DRY lambda = 0'
%%% in Figure 5 in B&K.
figure;
plot(S, Depth, 'r', 'LineWidth', 2)
grid on
ax = gca;
ax.YDir = 'reverse';
ax.FontSize = 14;
ax.XLim = [0 10];
ax.YLim = [0 5];
title('Fig. 2 from Dombard & McKinnon (2006)')
xlabel('Differential Stress Magnitude (MPa)', 'FontSize', 16)
ylabel('Depth (km)', 'FontSize', 16)