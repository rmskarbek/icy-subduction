function [Shear, mu, sigma_1] = IceShearResistance(Out, mu_flag)

%%% This function calculates the shear stress that resists sliding on an existing 
%%% fault surface with a variable dip angle. For IcySubduction.m, we assume that the 
%%% fault surface exists where the upper surface of the subducting slab is in contact 
%%% with the adjacent convective layer of the ice shell.

%%% The calculation follows Mueller & Phillips (1991). Along the fault surface, 
%%% brittle and ductile shear stresses are determined at a particular depth, and the 
%%% lesser of the two is assumed to represent the local shear stress resisting 
%%% sliding motion.

%%% The stress is integrated along the fault plane to obtain the resisting force 
%%% oriented along the fault. 
 
%%% The stress is also integrated with depth in the ice slab to obtain the resisting 
%%% force oriented in the direction of the first principal stress (in this case,
%%% horizontally).

%%% Finally, for comparison the resisting force is also calculated following Howell &
%%% Pappalardo (2019), who assume a constant dip angle for their friction force
%%% calculation, even though they use a variable dip angle in their subduction model.
%%% Also, H2019 assume a constant lithostatic stress in the slab instead of using a 
%%% depth dependent stress.

%%%-------------------------------------------------------------------------------%%%
%%% INPUT
%%%-------------------------------------------------------------------------------%%%
% Out      - A structure of simulation output generated from IcySubduction.m.

% mu_flag  - Determines how to handle the friction coefficient. 
%            Set mu_flag = 'temperature' to use temperature dependent values from 
%            Persson (2015). Set mu_flag = 'constant' to use a constant value.


%%%-------------------------------------------------------------------------------%%%
%%% OUTPUT. All of these variables are stored in a structure called 'Shear'.
%%%-------------------------------------------------------------------------------%%%
% S_Shear       - The shear resistance along the fault due friction and/creep.

% Tau_F         - The frictional shear resistance along the fault.

% Sigma_Diff    - The resistance along the fault due to creep.

% F_fault       - The integrated force along the fault, calculated from Shear.

% F_horizontal  - The integrated force along the fault in the x-direction, calculated 
%                 from Shear.

% F_frict       - The horizontal force from Tau_F, calculated following H2019 by 
%                 assuming a constant dip angle.

% Xi            - Along fault coordinates corresponding to the input depth values.


%%%-------------------------------------------------------------------------------%%%
%%% Geometry.
%%%-------------------------------------------------------------------------------%%%
p = Out.p;
H = p.Geometry.SlabThick;                       % thickness of conductive slab [m]

%%% Arc length of the upper slab surface, will be used as the along-fault distance.
Z_Top = Out.Slab.Top;
dydx = gradient(imag(Z_Top))./gradient(real(Z_Top));               % [-]
ds = sqrt(1 + dydx.^2);
Xi = cumtrapz(real(Z_Top), ds);                 % Along-fault distance [m].

%%% Use this for Xi to compare with H2019 equation, and set constant dip angle.
% Xi = Depth/sin(Theta);

%%% Depth of the upper slab surface.
Depth = Out.Slab.Depth(1,:)';
Theta = Out.Slab.Dip;                           % slab dip angle [radian]

%%% Temperature at the upper slab surface.
TempK = Out.Temperature(1,:)';

%%% Lithostatic stress at the upper slab surface. Here we assume that sigma_3 = sigma_L.
sigma_L = 1e-6*Out.VerticalStress(1,:)';

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
S_Shear = min(S_Shear,[],2);                                % shear resistance [MPa]

%%% Where theta = 0, change the shear resistance to zero.
i_0 = Theta == 0;
S_Shear(i_0) = 0;

%%% The plate interface does not exist at depths greater than the slab thickness, so 
%%% set S_Shear = nan there.
[~, i_H] = min(abs(Depth - H));
S_Shear(i_H+1:end) = nan;

%%% Calculate the force per unit distance directed along the fault plane by 
%%% integrating the shear resistance along the fault plane. Only integrate through 
%%% the conductive section of the ice shell.
F_Fault = 1e6*cumtrapz(Xi, S_Shear);                              % [N / m]

%%% Calculate the force per unit distance in the sigma_1 (horizontal) direction by
%%% integrating the additional horizontal stress that is needed to overcome the shear
%%% resistance. This integral is done along the vertical coordinate, i.e. with depth. 
%%% Only integrate through the conductive section of the ice shell.
S_Frict = (2./sin(2*Theta)).*S_Shear;
F_Horizontal = [1e6*cumtrapz(Depth(2:end), S_Frict(2:end)); nan];        % [N / m]

%%% H2019 frictional resisting force, equation (B12) in that paper. Here we use a 
%%% depth dependent lithostatic stress instead of a constant value like H2019 used.
mu_k = 0.37;
theta = pi*30.6/180;
sigma_3 = sigma_L;
sigma_Frict = mu_k*sigma_3/(cos(theta)*sin(theta) - mu_k*sin(theta)^2);
F_Frict = 1e6*cumtrapz(Depth(1:i_H), sigma_Frict(1:i_H));

%%% This is equation (B12) from H2019, where they use sigma_3 = g*rho_ice*H/2.
% F_frict = g*mu_k*rho_ice*H^2/(2*(cos(theta)*sin(theta) - mu_k*sin(theta)^2));

%%%-------------------------------------------------------------------------------%%%
%%% Output
%%%-------------------------------------------------------------------------------%%%
%%% Pad the force arrays with nan so they can easily be plotted against the arc 
%%% length, ect.
F_Frict = [F_Frict; nan(numel(Xi) - i_H,1)];

Shear.S_Shear = S_Shear;
Shear.S_Frict = S_Frict;
Shear.Tau_F = Tau_F;
Shear.Sigma_Diff = Sigma_Diff;
Shear.F_Fault = F_Fault;
Shear.F_Horizontal = F_Horizontal;
Shear.F_Frict = F_Frict;
Shear.Xi = Xi;
Shear.mu_flag = mu_flag;

end