function [Depth, DipAngle, Resistance, S_30] = FailureEnvelopes2

%%% Some constants from Dombard & McKinnon (2006).
D = 5;                                                      % shell thickness [km]
g = 1.31;                                                   % Europa gravity [m / s^2]
R = 8.314;                                                  % gas constant [J / (mol K]
rho_ice = 950;                                              % ice density [kg/m^3]
Tc = 273.15;                                                % water melt temp 1 atm [K]

%%% Dip angle of subduction fault plane [radians]. See Fig. 1 in Mueller & Phillips.
DipAngle = (pi/180)*(10:0.1:45)';
[~, i_30] = min(abs((180/pi)*DipAngle - 30));

%%% Set up a figure for plotting the results.
font = 'DejaVu Serif';
figure('DefaultTextFontName', font, 'DefaultAxesFontName', font,...
    'DefaultAxesFontSize',  14)
set(gcf, 'Color', 'w')

%%% Force per unit distance to cause sliding on the fault plane.
Resistance = nan(size(DipAngle));

%%% Flags for how to handle the friction coefficient and line styles for each case.
Mu_Flags = {'temperature'; 'constant'};
LineStyle = {'-'; '-.'};
for j = 1:numel(Mu_Flags)
    for i = 1:numel(DipAngle)
        [Xi, Depth, TempK, S, F] = ShearResistance(DipAngle(i), Mu_Flags{j});
        Resistance(i) = F;
    
        if i == i_30
        %%% Plot the shear resistance against the depth on the fault plane.
            subplot(1, 2, 1);
            hold on
            S_30 = S;
            plot(S, Depth, LineStyle{j}, 'LineWidth', 2)
            hold off
        end
    end
end

%%% Labels and such for the shear resistance plot.
box on
grid on
ax = gca;
ax.YDir = 'reverse';
ax.FontSize = 14;
title('Europa Shear Resistance, 30 deg dip')
xlabel('Shear Resistance (MPa)', 'FontSize', 16)
ylabel('Depth (km)', 'FontSize', 16)
legend('\mu = 0.55, Howell & Pappalardo (2019)', 'Temperature dependent');

%%% Plot the temperature profile next to the shear resistance profile.
subplot(1, 2, 2)
plot(TempK, Depth, 'k-', 'LineWidth', 2)
grid on
ax = gca;
ax.YDir = 'reverse';
ax.FontSize = 14;
% title('Europa Shear Resistance, 30 deg dip')
xlabel('Temperature (k)', 'FontSize', 16)
ylabel('Depth (km)', 'FontSize', 16)

%%% Plot the resistance against the dip angle.
% figure
% plot((180/pi)*DipAngle, 1e-9*Resistance, 'LineWidth', 2)
% grid on
% ax = gca;
% ax.FontSize = 14;
% ax.XLim = [0 50];
% % ax.YLim = [1e13 16e13];
% title('Europa, Shear Resitance Force')
% xlabel('Fault Plane Dip (deg)', 'FontSize', 16)
% ylabel('Resistance (GN/m)', 'FontSize', 16)

%%% Output dip angle in degrees.
DipAngle = (180/pi)*DipAngle;

function [Xi, Depth, TempK, S, F] = ShearResistance(theta, mu_flag)
    %%% Fault length [km].
    L = D/sin(theta);
    
    %%% Along-dip distance [km].
    Xi = linspace(0, L, 1000)';
    
    %%% Depth along the fault plane [km].
    Depth = Xi*sin(theta);

    %%% Temperature profile from Dombard & McKinnon (2006).
    T_s = 80;                                                   % Surface temperature [K]
    dTdz = 18;                                                  % thermal gradient [K/km]
    TempK = T_s + dTdz*Depth;                                   % [K]
    TempC = TempK - Tc;                                         % [C]
    
    %%% Ice density [kg/m^3] as a function of temperature. From Fukusako (1990),
    %%% Thermophysical Properties of Ice, Snow, and Sea Ice. There may be a more up to date
    %%% reference for ice density. However, the temperature dependence does not have much
    %%% affect relative to using a constant density.
    i_rho = TempC < -140;
    rho = 917*(1 - 1.17e-4*TempC);                              % eq. (4) in Fukusako
    rho(i_rho) = 930*(1 - 1.54e-5*TempC(i_rho));                % eq. (5) in Fukusako

    %%% Lithostatic stress along the fault plane. Here we assume that sigma_3 = sigma_L.
    sigma_L = 1e-3*g*cumtrapz(Depth, rho);                      % [MPa]    

    %%% Coefficient of friction from Persson (2015). 
    %%% 9.6.2024 - TempK_P is close to, but not exactly the same as TempK. So that needs to be
    %%% addressed.
    switch mu_flag
        case 'constant'
            mu = 0.55;                                      % Howell & Pappalardo, static
            % mu = 0.37;                                      % Howell & Pappalardo, kinetic

        case 'temperature'
            [TempK_P, mu] = IceFriction(TempK);

    end    
    C = 0;                                                      % Cohesion [MPa]

    %%% Compute sigma_1 in the brittle regime using the friction law, from Mueller & Phillips 
    %%%, eq (5). This is not a failure envelope.
    sigma_1 = (2*C + (sin(2*theta) + mu*(cos(2*theta) + 1)).*sigma_L)...
        ./(sin(2*theta) + mu*(cos(2*theta) - 1));
    
    %%% Frictional shear resistance on the fault plane in terms of principal stresses.
    %%% Also not a failure envelope calculation.
    Tau_F = 0.5*(sigma_1 - sigma_L)*sin(2*theta);
    
    %%% Ice flow laws, paramters for Regimes A, B, and C are from from Dombard & McKinnon 
    %%% (2001), Table 1. These regimes represent grainsize insensitive disclocation creep 
    %%% mechanisms. See Durham & Stern (2001) for references. See text in D & M (2006) in the
    %%% paragraph after eq (1b). Symbols here are after Dombard & McKinnon (2001).
    e = 1e-15;                                              % strain rate [1/s]
    
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
    
    %%% Find the shear resistance.
    S = [Tau_F, Tau_A, Tau_B, Tau_C, Tau_G];
    S = min(S,[],2);
    
    %%% Minimum force (pre unit distance) to promote incipient convergence is found by
    %%% integrating the shear resistance along the fault plane.
    F = 1e9*trapz(Xi, S);                                   % [N / m]

end
end