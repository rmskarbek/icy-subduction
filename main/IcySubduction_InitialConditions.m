function [Vars_init, sigma_lith] = IcySubduction_InitialConditions(p)

%%% Initial conditions for IcySubduction.m. This function is called by 
%%% IcySubduction_Parameters.m, and follows the procedure detailed in Howell2019, 
%%% Appendix A.1, to generate a self-consisent porosity profile.

%%%-------------------------------------------------------------------------------%%%
%%% Constants.
%%%-------------------------------------------------------------------------------%%%
    % spy = 365.25*24*3600;                   % seconds per year
    g = p.Constants.Gravity;                % Europa gravity [m/s^2]
    R = p.Constants.GasConst;               % gas constant [J/(mol K)]
    
%%%-------------------------------------------------------------------------------%%%
%%% Geometry
%%%-------------------------------------------------------------------------------%%%
    H = p.Geometry.SlabThick;               % thickness of conductive slab [m]
    N = p.Numerical.GridPoints;             % Number of grid points

%%%-------------------------------------------------------------------------------%%%
%%% Temperature
%%%-------------------------------------------------------------------------------%%%
    T_b = p.Temperature.BasalTemp;          % basal temp [K] 
    T_s = p.Temperature.SurfaceTemp;        % surface temp [K], after Nimmo et al. (2003)

%%%-------------------------------------------------------------------------------%%%
%%% Viscosity
%%%-------------------------------------------------------------------------------%%%
    Q = p.Viscosity.ActivEnergy;            % viscosity activation energy [J/mol]
    eta_b = p.Viscosity.BasalVisc;          % viscosity in convecting ice T = T_b [Pa*s]
    eta_max = p.Viscosity.MaxVisc;          % maximum viscosity from Howell2019 [Pa*s]

%%%-------------------------------------------------------------------------------%%%
%%% Porosity and Density
%%%-------------------------------------------------------------------------------%%%
    t_por = p.Porosity.CompTime;            % initial porosity compaction time [s]
    phi_0 = p.Porosity.PorosityRef;         % porosity reference
    f_slab = p.Porosity.SaltSlab;           % salt content in slab
    rho_salt = p.Porosity.SaltDensity;      % salt density [kg/m^3]
    method = p.Porosity.DensityMethod;      % method for computing ice density

%%%-------------------------------------------------------------------------------%%%
%%% Compute the initial temperature and porosity profiles.
%%%-------------------------------------------------------------------------------%%%
%%% Depth referenced to the top of the column [m].
    Depth_init = linspace(0, H, N)';

%%% Initial temperature profile in the subducting slab. Eq. (3) in Johnson2017, [K].
    T_init = T_s*(T_b/T_s).^(Depth_init/H);

%%% Initial viscosity profile in the subducting slab. Eq. (5) in Johnson2017, [Pa s].
    eta_init = eta_b*exp((Q/R)*(1./T_init - 1/T_b));
    eta_init(eta_init > eta_max) = eta_max;

%%% Initial lithostatic stress [Pa].
    rho_ice = IceDensity(T_init, method);
    rho = (1 - f_slab)*rho_ice + f_slab*rho_salt;
    sigma_lith = g*rho.*Depth_init;    

%%% Iterate through the porosity calculation until the lithostatic stress stops 
%%% changing by some amount.
    tol_stress = 1e-10;             % [Pa]
    q = 0;
    t = 1;
    while t > tol_stress
        q = q + 1;
        sigma_0 = sigma_lith;

%%% Initial porosity profile in the subducting slab. Eq. (4) in Johnson2017
        phi_init = phi_0*exp(-t_por*sigma_lith./eta_init);
        rho = ((1 - f_slab)*rho_ice + f_slab*rho_salt).*(1 - phi_init);
        sigma_lith = cumtrapz(Depth_init, g*rho);

        t = max(abs(sigma_lith - sigma_0));
    end

%%% Construct the intial conditions vector.
    Vars_init = [T_init; phi_init];

end