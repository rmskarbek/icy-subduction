function [Buoyancy, delta_rho_ice, F_Buoy_salty]...
    = BuoyancySalt(Out, Bend, Buoyancy, Shear)

%%% This function takes simulation output and determines how much denser the ice
%%% would need to be, in order to offset the resisting forces F_shear and F_bend.


%%% Find the arc length location where the slab is completely subsumed.
[~, ~, ~, k_subsumed] = SlabSubsumption(Out);

%%% Total resistive force.
i_f = isfinite(Shear.F_Horizontal);
S_f = Shear.F_Horizontal(i_f);
F_Resist = Bend.F_Fail(k_subsumed) + S_f(end);

%%% Run a bisection method search to find the increased ice grain density where the 
%%% buoyancy force is equal to F_bend + F_shear. Assume that this occurs somewhere in
%%% the interval delta_rho/rho_ice_0 = [a b]. Where rho_ice_0 is the density of pure
%%% ice at 0 degrees Celsius.
rho_ice_0 = IceDensity(273.15, 'density');
a = 0;
b = 3;

%%% Search tolerance [kg/m^3].
tol = 1;

%%% Maximum number of steps.
r_max = 100;
r = 1;
while r < r_max

    %%% Stop the loop if the half the length of the current interval is less
    %%% than the tolerance.
    if ((b-a)/2)*rho_ice_0 < tol
        break
    end

    %%% Check the buoyance force for the increased ice grain density.
    delta_hat = (a + b)/2;
    delta_rho_ice = delta_hat*rho_ice_0;
    Buoyancy_salty = BuoyancyForce(Out, delta_rho_ice);
    F_Buoy_salty = Buoyancy_salty.F_Bouy(k_subsumed);

    %%% Update the search interval.       
    if F_Buoy_salty < F_Resist
        a = delta_hat;
    else
        b = delta_hat;
    end

    %%% Increment the step counter.
    r = r + 1;
end

%%% Store the computed delta_rho_ice.
Buoyancy.Delta_Rho = delta_rho_ice;
Buoyancy.F_Buoy_Salty = F_Buoy_salty;