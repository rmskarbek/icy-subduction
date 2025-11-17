function p = IcySubduction_Parameters
%%% This function generates a structure of parameter values and constants that are 
%%% used in IcySubduction.m. Having a separate code for defining the parameters helps
%%% in maintaining records of simulation data, and in running large batches of 
%%% simulations.

%%% This function calls IcySubduction_InitialConditions.m to generate the initial 
%%% conditions for IcySubduction.m

%%%-------------------------------------------------------------------------------%%%
%%% Constants.
%%%-------------------------------------------------------------------------------%%%
spy = 365.25*24*3600;               % seconds per year
g = 1.31;                           % Europa gravity [m/s^2]
R = 8.314;                          % gas constant [J/(mol K)]
% rho_ice0 = 917;                     % reference ice density at T = 0 C[kg/m^3]

%%%-------------------------------------------------------------------------------%%%
%%% Geometry
%%%-------------------------------------------------------------------------------%%%
% H = 2.5e3;                          % thickness of conductive slab [m]
H = 5e3;                            % thickness of conductive slab [m]
H_shell = 10e3;                     % total thickness of ice shell [m]
GeoFlag = 'Buffett';                % a flag that determines the slab geometry.
% GeoFlag = 'Johnson';                % a flag that determines the slab geometry.
R_min = 2*H;                        % minimum radius of plate curvature [m]
v_plate = 0.004;                     % plate convergence rate [m/year]

%%%-------------------------------------------------------------------------------%%%
%%% Numerical
%%%-------------------------------------------------------------------------------%%%
N = 200;                            % Number of grid points.
dz = H/(N-1);                       % grid spacing [m]

%%% Flag for how to handle variable thermal conductivity in the finite difference 
%%% scheme.
% coeff = 'linear';                 % use a linear mean, see Zikanov
coeff = 'harmonic';                 % use a harmonic mean

%%% Compute the simulation run time so that the simulation ends when the bottom edge 
%%% of the slab reaches the bottom of the ice shell.
type = 'geometry';
[~, ~, ~, ~, ~, ~, s_end] = Buffett2006(type, H, H_shell, R_min);
RunTime = spy*s_end/v_plate;

%%%-------------------------------------------------------------------------------%%%
%%% Temperature
%%%-------------------------------------------------------------------------------%%%
T_b = 260;                          % basal temp [K] 
T_s = 100;                          % surface temp [K]after Nimmo et al. (2003)

%%%-------------------------------------------------------------------------------%%%
%%% Viscosity
%%%-------------------------------------------------------------------------------%%%
Q = 50000;                          % viscosity activation energy [J/mol]
eta_b = 1e13;                       % reference viscosity at base of ice shell [Pa*s]
eta_max = 1e23;                     % maximum viscosity from Howell2019 [Pa*s]

%%%-------------------------------------------------------------------------------%%%
%%% Porosity and Density
%%%-------------------------------------------------------------------------------%%%
t_por = 65e6;                       % initial porosity compaction time [years]
phi_0 = 0.1;                        % porosity reference
% phi_0 = 0;                          % porosity reference

%%% Salt content in the slab, non-subducting conductive layer, and the convective 
%%% layer.
rho_salt = 1444;                    % natron density [kg/m^3]

% f_slab = 0;                         % salt content in slab
f_conduct = 0;                      % salt content in convective layer
f_convect = 0;                      % salt content in non-subducting conductive layer

f_slab = 0.1;                         % salt content in slab
% f_conduct = 0.01;                   % salt content in convective layer
% f_convect = 0.15;                   % salt content in non-subducting conductive layer


%%% Method for computing ice density; see IceDensity.m
% method = 'constant';                % use a constant value, defined in IceDensity.m
method = 'density';                 % temperature dependent, from ice density data.
% method = 'thermal';                 % temperature dependent, from thermal expansion.

%%%-------------------------------------------------------------------------------%%%
%%% Parameters Structure
%%%-------------------------------------------------------------------------------%%%
p = struct('Constants', [], 'Numerical', [], 'Geometry', [], 'Temperature', [],...
    'Viscosity', [], 'Porosity', []);

p.Constants.Gravity = g;
p.Constants.GasConst = R;
% p.Constants.IceDensity = rho_ice0;

p.Geometry.SlabThick = H;
p.Geometry.ShellThick = H_shell;
p.Geometry.CurveRadius = R_min;
p.Geometry.GeoFlag = GeoFlag;
p.Geometry.PlateRate = v_plate/spy;     % convert to [m/s]

p.Numerical.GridPoints = N;
p.Numerical.GridSpacing = dz;
p.Numerical.RunTime = RunTime;
p.Numerical.VarCoeff = coeff;

p.Temperature.SurfaceTemp = T_s;
p.Temperature.BasalTemp = T_b;

p.Viscosity.ActivEnergy = Q;
p.Viscosity.BasalVisc = eta_b;
p.Viscosity.MaxVisc = eta_max;

p.Porosity.CompTime = spy*t_por;        % convert to [s]
p.Porosity.PorosityRef = phi_0;

p.Porosity.SaltDensity = rho_salt;
p.Porosity.SaltSlab = f_slab;
p.Porosity.SaltConduct = f_conduct;
p.Porosity.SaltConvect = f_convect;

p.Porosity.DensityMethod = method;

%%%-------------------------------------------------------------------------------%%%
%%% Initial conditions
%%%-------------------------------------------------------------------------------%%%
Vars_Init = IcySubduction_InitialConditions(p);
p.Numerical.InitialCond = Vars_Init;

end