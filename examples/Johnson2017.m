%%% Parameters for the "fiducial" simulation from Johnson et al. (2017).
H_slab = 6e3;
Radius = 2*H_slab;
Dip = 45;
RunTimeYRS = 550e3;
% method = 'constant';
method = 'thermal';
eta_b = 1e15;

%%% Generate parameters structure and update it with the Johnson parameters.
p = IcySubduction_Parameters;

p.Geometry.SlabThick = H_slab;
p.Geometry.CurveRadius = Radius;
p.Geometry.DipAngle = Dip;
p.Geometry.GeoFlag = 'CircularArc';
% p.Viscosity.BasalVisc = eta_b;
p.Viscosity.MaxVisc = inf;
p.Porosity.PorosityRef = 0.1;
p.Porosity.DensityMethod = method;
% p.Porosity.SaltContent = 0.02;

spy = 365.25*24*3600;   % seconds per year
p.Numerical.RunTime = spy*RunTimeYRS;

%%% Update the initial conditions. Johnson2017 did not use the iterative method that
%%% Howell2019 used.
Vars_Init = IcySubduction_InitialConditions(p);
p.Numerical.InitialCond = Vars_Init;

%%% Run the simulations and plot the output after Figures 3 and 4 in Johnson2017.
Out = IcySubduction(p);
Bouyancy = BouyancyForce(Out);

Johnson2017_Plots1(Out, Bouyancy.Rho_Anom);