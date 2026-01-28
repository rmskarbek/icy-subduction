function Buoyancy = BuoyancyForce(Out)

%%% This function computes the density anomaly by comparing the densities in the slab
%%% with those outside of the slab at equivalent depths.

%%%-------------------------------------------------------------------------------%%%
%%% Geometry.
%%%-------------------------------------------------------------------------------%%%
p = Out.p;
g = p.Constants.Gravity;
H = p.Geometry.SlabThick;               % thickness of conductive slab [m]

%%%-------------------------------------------------------------------------------%%%
%%% Simulation Output.
%%%-------------------------------------------------------------------------------%%%
Time = Out.TimeYrs;
ArcLength = Out.Slab.ArcLength;
Theta = Out.Slab.Dip;
Depth = Out.Slab.Depth;
Rho = Out.Density;

%%%-------------------------------------------------------------------------------%%%
%%% Densities outside of the slab at equivalent depths.
%%%-------------------------------------------------------------------------------%%%
Rho_Shell = nan(size(Rho));             % Bulk density in the ice shell
f_conduct = p.Porosity.SaltConduct;     % salt content in non-subducting conductive layer
f_convect = p.Porosity.SaltConvect;     % salt content in convective layer

%%% Find depths that are less than the conductive layer thickenss.
i_Cond = Depth <= H;
i_Conv = i_Cond ~= 1;

%%% Ice grain density in the convecting layer.
T_b = p.Temperature.BasalTemp;          % basal temp [K] 
method = p.Porosity.DensityMethod;      % method for computing ice density
rho_ice_b = IceDensity(T_b, method);    % density of convecting ice.
rho_salt = p.Porosity.SaltDensity;      % salt density [kg / m^3]

%%% Bulk density in the convecting layer - is a constant.
Rho_Shell(i_Conv) = (1 - f_convect)*rho_ice_b + f_convect*rho_salt;

%%% Ice grain density in the non-subducting conductive layer.
Temp_init = Out.Temperature(:,1);
phi_init = Out.Porosity(:,1);
rho_ice_init = IceDensity(Temp_init, method);
rho_cond = ((1 - f_conduct)*rho_ice_init + f_conduct*rho_salt).*(1 - phi_init);

%%% Bulk density in the non-subducting conductive layer.
Rho_Shell(i_Cond) = interp1(Depth(:,1), rho_cond, Depth(i_Cond));

%%% At any location (s,z) in the column, the density anomaly is the difference
%%% between the slab density at (s,z), and the ice shell density at the same depth z.
%%% Then take the mean over the whole column.
Rho_Anom = mean(Rho - Rho_Shell)';

%%%-------------------------------------------------------------------------------%%%
%%% Buoyancy Force.
%%%-------------------------------------------------------------------------------%%%
%%% Incorrect equation (B1) from Howell2019, does not account for slab dip.
% F_Buoy = zeros(size(Time));
% F_Buoy(2:end) = g*Rho_Anom(2:end).*H.*diff(ArcLength);
% F_Buoy = cumsum(F_Buoy);                % [N / m]

%%% Correct equation from Appendix A in Li and Gurnis, does account for slab dip.
F_int = g*H*sin(Theta).*Rho_Anom;
F_Buoy = cumtrapz(ArcLength, F_int);

%%%-------------------------------------------------------------------------------%%%
%%% Output.
%%%-------------------------------------------------------------------------------%%%
Buoyancy.Rho_Shell = Rho_Shell;
Buoyancy.Rho_Anom = Rho_Anom;
Buoyancy.F_Bouy = F_Buoy;