function Shear = IceShearResistance(Out, mu_flag)

%%% This function calculates the shear stress that resists sliding on an existing fault
%%% surface with a variable dip angle. For IcySubduction.m, we assume that the fault 
%%% surface exists where the upper surface of the subducting slab is in contact with the
%%% adjacent convective layer of the ice shell; this is indicated by the red, dashed line
%%% in the figure produced by GeometryPlot.m.

%%% The calculation follows Mueller & Phillips (1991). Along the fault surface, brittle
%%% and ductile shear stresses are determined at a particular depth, and the lesser of the
%%% two is assumed to represent the local shear stress resisting sliding motion.

%%% The stress is integrated along the fault plane to obtain the resisting force oriented
%%% along the fault. The stress is also integrated with depth in the ice slab to obtain
%%% the resisting force oriented in the direction of the first principal stress (in this
%%% case, horizontally).

%%% Finally, for comparison the resisting force is also calculated following Howell &
%%% Pappalardo (2019), who assume a constant dip angle for their friction force
%%% calculation, even though they use a variable dip angle in their subduction model.
%%% Also, H2019 assume a constant lithostatic stress in the slab instead of using a depth
%%% dependent stress.

%%% 9.18.2024 - Input cohesion? Then can run a test that duplicates DombardMcKinnon2006.m

%%% 9.19.2024 - Shear resistance becomes negative if the dip angle is greater than 45
%%% degrees?

%%% 3.11.2024 - sigma_1 = inf, for theta = 0?

%%%------------------------------------------------------------------------------------%%%
%%% INPUT
%%%------------------------------------------------------------------------------------%%%
% mu_flag  - Determines how to handle the friction coefficient. Set mu_flag = 'temperature'
%            to use temperature dependent values from Persson (2015).
%            Set mu_flag = 'constant' to use a constant value.

%%%------------------------------------------------------------------------------------%%%
%%% OUTPUT. All of these variables are stored in a structure called 'Shear'.
%%%------------------------------------------------------------------------------------%%%
% S_Shear       - The shear resistance along the fault due friction and/creep.

% Tau_F         - The frictional shear resistance along the fault.

% Sigma_Diff    - The resistance alonf the fault due to creep.

% F_fault       - The integrated force along the fault, calculated from Shear.

% F_horizontal  - The integrated force along the fault in the x-direction, calculated from 
%                 Shear.

% F_frict       - The horizontal force from Tau_F, calculated following H2019 by assuming
%                 a constant dip angle.

% Xi            - Along fault coordinates corresponding to the input depth values.


%%%------------------------------------------------------------------------------------%%%
%%% Constants.
%%%------------------------------------------------------------------------------------%%%
p = Out.p;
g = p.Constants.Gravity;                        % Europa gravity [m/s^2]
% Tc = 273.15;                                    % water melt temp 1 atm [K]
rho_ice = p.Constants.IceDensity;               % reference ice density [kg/m^3]

%%%------------------------------------------------------------------------------------%%%
%%% Geometry.
%%%------------------------------------------------------------------------------------%%%
H = p.Geometry.SlabThick;                       % thickness of conductive slab [m]

%%% Arc length of the upper slab surface, will be used as the along-fault distance.
Z_Top = Out.Slab.Top;
dydx = gradient(imag(Z_Top))./gradient(real(Z_Top));               % [-]
ds = sqrt(1 + dydx.^2);
Xi = cumtrapz(real(Z_Top), ds);                 % Along-fault distance [m].
% Xi = Out.Slab.ArcLength;

%%% Use this for Xi to compare with H2019 equation, and set constant dip angle.
% Xi = Depth/sin(Theta);

%%% Depth of the upper slab surface.
Depth = Out.Slab.Depth(1,:)';
Theta = Out.Slab.Dip;                           % slab dip angle [radian]

% Z_Center = Out.CenterLine;
% dydx_Center = gradient(imag(Z_Center))./gradient(real(Z_Center));        % [-]
% Theta = atan(dydx_Center);                                               % [radian]
% Theta = (pi/180)*30.6;

%%% Temperature at the upper slab surface.
TempK = Out.Temperature(1,:)';

%%% Lithostatic stress at the upper slab surface. Here we assume that sigma_3 = sigma_L.
sigma_L = 1e-6*g*rho_ice*Depth;                 % [MPa]

%%% Coefficient of friction.
switch mu_flag
    case 'constant'
        % mu = 0.55;                            % Howell & Pappalardo, static
        mu = 0.37;                              % Howell & Pappalardo, kinetic

    case 'temperature'
    %%% Coefficient of friction from Persson (2015).        
        [TempK_P, mu] = IceFriction(TempK);
end

%%% Cohesion [MPa].
C = 0;
% C = 1;

%%% Compute sigma_1 in the brittle regime using the friction law.
sigma_1 = (2*C + (sin(2*Theta) + mu.*(cos(2*Theta) + 1)).*sigma_L)...
    ./(sin(2*Theta) + mu.*(cos(2*Theta) - 1));

%%% Frictional shear resistance on the fault/failure plane in terms of principal 
%%% stresses.
Tau_F = 0.5*(sigma_1 - sigma_L).*sin(2*Theta);
        
%%% Ductile differential stress from ice flow laws.
Sigma_Diff = IceFlowLaws(TempK);

%%% Find the shear resistance.
S_Shear = [Tau_F, Sigma_Diff];
S_Shear = min(S_Shear,[],2);                                    % shear resistance [MPa]

%%% Where theta = 0, change the shear resistance to nan.
i_nan = Theta == 0;
S_Shear(i_nan) = 0;

%%% Calculate the force per unit distance directed along the fault plane by integrating
%%% the shear resistance along the fault plane. Only integrate through the conductive
%%% section of the ice shell.
[~, i_H] = min(abs(Depth - H));
% F_fault = 1e6*trapz(Xi(1:i_cond), S_Shear(1:i_cond));                   % [N / m]
F_Fault = 1e6*cumtrapz(Xi(1:i_H), S_Shear(1:i_H));                % [N / m]

% F_fault = 1e6*trapz(Xi(1:i_cond), Tau_F(1:i_cond));                     % [N / m]

%%% Calculate the force per unit distance in the sigma_1 (horizontal) direction by
%%% integrating the shear resistance with depth. Only integrate through the conductive
%%% section of the ice shell.
% F_horizontal = 1e6*trapz(Depth(1:i_cond), S_Shear(1:i_cond));           % [N / m]
F_Horizontal = 1e6*cumtrapz(Depth(1:i_H), S_Shear(1:i_H));        % [N / m]

% F_horizontal = 1e6*trapz(Depth, Tau_F);

%%% H2019 frictional resisting force, equation (B12) in that paper. Here we use a depth
%%% dependent lithostatic stress instead of a constant value like H2019 used.
mu_k = 0.37;
theta = pi*30.6/180;
sigma_3 = sigma_L;
sigma_Frict = mu_k*sigma_3/(cos(theta)*sin(theta) - mu_k*sin(theta)^2);
% F_frict = 1e6*trapz(Depth(1:i_cond), sigma_frict(1:i_cond));
F_Frict = 1e6*cumtrapz(Depth(1:i_H), sigma_Frict(1:i_H));

% m * MPa = (MN / m) * 1e6 N / MN

%%% This is equation (B12) from H2019, where they use sigma_3 = g*rho_ice*H/2.
% F_frict = g*mu_k*rho_ice*H^2/(2*(cos(theta)*sin(theta) - mu_k*sin(theta)^2));

%%%------------------------------------------------------------------------------------%%%
%%% Output
%%%------------------------------------------------------------------------------------%%%
%%% Pad the force arrays with nan so they can easily be plotted against the arc length,
%%% ect.
F_Fault = [F_Fault; nan(numel(Xi) - i_H,1)];
F_Horizontal = [F_Horizontal; nan(numel(Xi) - i_H,1)];
F_Frict = [F_Frict; nan(numel(Xi) - i_H,1)];

Shear.S_Shear = S_Shear;
Shear.Tau_F = Tau_F;
Shear.Sigma_Diff = Sigma_Diff;
Shear.F_Fault = F_Fault;
Shear.F_Horizontal = F_Horizontal;
Shear.F_Frict = F_Frict;
Shear.Xi = Xi;
Shear.mu_flag = mu_flag;

end

%%% Ice density [kg/m^3] as a function of temperature. From Fukusako (1990),
%%% Thermophysical Properties of Ice, Snow, and Sea Ice. There may be a more up to date
%%% reference for ice density. However, the temperature dependence does not have much
%%% affect relative to using a constant density.
% TempC = TempK - Tc;                                         % [C]
% i_rho = TempC < -140;
% rho = 917*(1 - 1.17e-4*TempC);                              % eq. (4) in Fukusako
% rho(i_rho) = 930*(1 - 1.54e-5*TempC(i_rho));                % eq. (5) in Fukusako
% rho = rho_ice*ones(numel(Depth),1);