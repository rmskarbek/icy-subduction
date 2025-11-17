function Buoyancy = BuoyancyForce(Out)

%%% This function computes the density anomaly by comparing the densities in the slab
%%% with those outside of the slab at equivalent depths.

%%%-------------------------------------------------------------------------------%%%
%%% Geometry
%%%-------------------------------------------------------------------------------%%%
p = Out.p;
g = p.Constants.Gravity;
H = p.Geometry.SlabThick;               % thickness of conductive slab [m]

%%%-------------------------------------------------------------------------------%%%
%%% Simulation Output.
%%%-------------------------------------------------------------------------------%%%
Time = Out.TimeYrs;
ArcLength = Out.Slab.ArcLength;
Depth = Out.Slab.Depth;
Rho = Out.Density;

%%%-------------------------------------------------------------------------------%%%
%%% Densities outside of the slab at equivalent depths.
%%%-------------------------------------------------------------------------------%%%
Rho_Shell = nan(size(Rho));
f_conduct = p.Porosity.SaltConduct;     % salt content in non-subducting conductive layer
f_convect = p.Porosity.SaltConvect;     % salt content in convective layer

%%% Find depths that are less than the conductive layer thickenss.
i_Cond = Depth <= H;
i_Conv = i_Cond ~= 1;

%%% Density in the convecting layer.
T_b = p.Temperature.BasalTemp;          % basal temp [K] 
method = p.Porosity.DensityMethod;      % method for computing ice density
rho_ice_b = IceDensity(T_b, method);    % density of convecting ice.
rho_salt = p.Porosity.SaltDensity;      % salt density [kg/m^3]

% Rho_Shell(i_Conv) = rho_ice_b;
Rho_Shell(i_Conv) = (1 - f_convect)*rho_ice_b + f_convect*rho_salt;

%%% Bulk density in the non-subducting conductive layer.
Temp_init = Out.Temperature(:,1);
phi_init = Out.Porosity(:,1);
rho_ice_init = IceDensity(Temp_init, method);
% rho_init = rho_ice_init.*(1 - phi_init);
rho_cond = ((1 - f_conduct)*rho_ice_init + f_conduct*rho_salt).*(1 - phi_init);
Rho_Shell(i_Cond) = interp1(Depth(:,1), rho_cond, Depth(i_Cond));

%%% Bulk density in the non-subducting conductive layer.
% Rho_Shell(i_Cond) = interp1(Depth(:,1), Rho_cond(:,1), Depth(i_Cond));

%%% Average density anomaly, after Johnson2017;
Rho_Anom = mean(Rho - Rho_Shell)';

F_Buoy = zeros(size(Time));
F_Buoy(2:end) = g*Rho_Anom(2:end).*H.*diff(ArcLength);
F_Buoy = cumsum(F_Buoy);

%%%-------------------------------------------------------------------------------%%%
%%% Output
%%%-------------------------------------------------------------------------------%%%
Buoyancy.Rho_Shell = Rho_Shell;
Buoyancy.Rho_Anom = Rho_Anom;
Buoyancy.F_Bouy = F_Buoy;
