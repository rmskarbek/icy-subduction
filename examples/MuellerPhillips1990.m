function MuellerPhillips1990
%%% This function reproduces Figures 2 and 3 in Mueller & Phillips (1990), On the
%%% Initiation of Subduction. The calculation consists of computing the force required to
%%% initiate motion on an existing fault plane with a given dip angle. This serves as an
%%% estimate of the minimum force necessary to promote trench formation, because other
%%% resisting forces (e.g. plate bending) are ignored.

%%% Some constants.
D = 75;                                                     % slab thickness [km]
g = 9.81;                                                   % Earth gravity                   
R = 8.314;                                                  % gas constant [J / (mol K
rho = 3000;                                                 % lithosphere density [kg/m^3]
rho_w = 1000;                                               % water density [kg/m^3]
% Tc = 273.15;                                                % water melt temp 1 atm [K]

%%% Dip angle of subduction fault plane [radians]. See Fig. 1 in Mueller & Phillips.
DipAngle = (pi/180)*(10:0.1:45)';
[~, i_30] = min(abs((180/pi)*DipAngle - 30));

%%% Force per unit distance to cause sliding on the fault plane.
Resistance = nan(size(DipAngle));

for i = 1:numel(DipAngle)
    [Xi, Depth, Tau_F_0, Tau_F_200, Tau_O, S, F] = ShearResistance(DipAngle(i));
    Resistance(i) = F;

    if i == i_30
    %%% Plot the shear resistance against the depth on the fault plane to reproduce Fig.2
    %%% in Muellet & Phillips for a 100 My old slab.
        figure;
        plot(S, Depth, 'r', 'LineWidth', 2)
        grid on
        ax = gca;
        ax.YDir = 'reverse';
        ax.FontSize = 14;
        title('Fig. 2 from Mueller & Phillips (1991). 100 My, 30 deg dip')
        xlabel('Shear Resistance (MPa)', 'FontSize', 16)
        ylabel('Depth (km)', 'FontSize', 16)
    end
end

%%% Plot the resistance against the dip angle to reproduce Fig. 3 in Mueller & Phillips
%%% for a 100 My old slab.
figure
plot((180/pi)*DipAngle, Resistance, 'LineWidth', 2)
grid on
ax = gca;
ax.FontSize = 14;
ax.XLim = [10 50];
ax.YLim = [1e13 16e13];
title('Fig. 3 from Mueller & Phillips (1991). 100 My')
xlabel('Fault Plane Dip (deg)', 'FontSize', 16)
ylabel('Resistance (N/m)', 'FontSize', 16)

function [Xi, Depth, Tau_F_0, Tau_F_200, Tau_O, S, F] = ShearResistance(theta)
    %%% Fault length [km].
    L = D/sin(theta);
    
    %%% Along-dip distance [km].
    Xi = linspace(0, L, 1000)';
    
    %%% Depth along the fault plane [km].
    Depth = Xi*sin(theta);
    
    %%% Lithostatic stress along the fault plane. Here we assume that sigma_3 = sigma_L.
    sigma_L = 1e-3*g*rho*Depth;                 % [MPa]
    
    %%% Hydrostatic pore fluid pressure along the fault plane.
    P_f = 1e-3*g*rho_w*Depth;
    
    %%% Compute sigma_1 in the brittle regime using the friction law, from Mueller & Phillips 
    %%%, eq (5). This is not a failure envelope.
    mu_0 = 0.85;                                % friction below 200 MPa normal stress
    mu_200 = 0.6;                               % friction above 200 MPa normal stress
    C_0 = 0;                                    % cohesion below 200 MPa normal stress [MPa]
    C_200 = 50;                                 % cohesion above 200 MPa normal stress [MPa]
    
    %%% Sigma_1 below 200 MPa normal stress.
    sigma_1_0 = (2*(C_0 - mu_0*P_f) + (sin(2*theta) + mu_0*(cos(2*theta) + 1)).*sigma_L)...
        ./(sin(2*theta) + mu_0*(cos(2*theta) - 1));
    
    %%% Sigma_1 above 200 MPa normal stress.
    sigma_1_200 = (2*(C_200 - mu_200*P_f) + (sin(2*theta) + mu_200*(cos(2*theta) + 1)).*sigma_L)...
        ./(sin(2*theta) + mu_200*(cos(2*theta) - 1));
    
    %%% Frictional shear resistance on the fault plane in terms of principal stresses.
    %%% Also not a failure envelope calculation.
    Tau_F_0 = 0.5*(sigma_1_0 - sigma_L)*sin(2*theta);
    Tau_F_200 = 0.5*(sigma_1_200 - sigma_L)*sin(2*theta);
    
    %%% Flow Laws.
    %%% Thermal parameters for half-space cooling from Table 1 in Mueller and Phillips.
    T_s = 270;                                  % surface temperature [K]
    T_a = 1650;                                 % asthenosphere temperature [K]
    kappa = 1.2e-6;                             % thermal diffusivity [m^2 / s]
    
    %%% Temperature profile. Half-space cooling model, eq (4.113) in Turcotte and Schubert.
    time = 100e6;                                               % age of slab [Myr]
    spy = 365.25*24*3600;
    TempK = T_a + (T_s - T_a)*erfc(1e3*Depth/(2*sqrt(kappa*spy*time)));
    
    %%% Wet olivine flow law. Parameters from Table 2 in Mueller & Phillips. In the text they
    %%% refer to this as "wet dunite". Symbols are as in Mueller & Phillips.
    H = 0.444e6;                                            % [J / mol]
    A = 7e4;                                                % [1 / (s MPa^3)]
    n = 3.4;
    e = 1e-15;                                              % strain rate [1/s]
    Tau_O = ((e/A)*exp(H./(R*TempK))).^(1/n);               % [MPa]
    
    %%% Find the shear resistance. Again, this is not a failure envelope!
    S = [Tau_F_0, Tau_F_200, Tau_O];
    S = min(S,[],2);        
    
    %%% Minimum force (pre unit distance) to promote incipient convergence is found by
    %%% integrating the shear resistance along the fault plane.
    F = 1e9*trapz(Xi, S);                                   % [N / m]

end
end