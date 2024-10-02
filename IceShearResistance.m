function [Depth, TempK, S, F, varargout] = IceShearResistance(Depth, TempK, calc, mu_flag,...
    varargin)
%%% This function calculates the shear resistance on an existing fault, or computes a
%%% failure envelope by assuming frictional failure takes place on optimally oriented
%%% planes.

%%% 9.18.2024 - Input cohesion? Then can run a test that duplicate DombardMcKinnon2006.m

%%% 9.19.2024 - Shear resistance becomes negative if the dip angle is greater than 45
%%% degrees?

%%% INPUT
% Temperature
% Depth
% type of calculation
% fault dip angle if assuming a preexisting fault
% flag for handling the friction coefficient

%%% OUTPUT
% Shear resistance or failure envelope
% Integrated force along fault coordinates if assuming a preexisting fault, or with depth
% for a failure envelope.

%%% Some constants.
    g = 1.31;                                               % Europa gravity [m / s^2]    
    rho_ice = 950;                                          % ice density [kg/m^3]
    Tc = 273.15;                                            % water melt temp 1 atm [K]
    TempC = TempK - Tc;                                     % [C]
    
%%% The shell thickness should be equal to the maximum input depth.
    % D = max(Depth);                                        % shell thickness [km]

%%% Handle the calculation type.    
    switch calc   
        case 'fault'
        %%% Calculation for a pre-existing fault. Check that theta has been input.
            if isempty(varargin)
                fprintf('For calc = ''fault'', input must be:\nIceShearResistance(Depth, TempK, calc, mu_flag, theta).\n');
                return
            end

        %%% Fault length [km].
            theta = (pi/180)*varargin{1};
            % L = D/sin(theta);
        
        %%% Along-fault distance [km].
            % Xi = linspace(0, L, numel(Depth))';
            Xi = Depth/sin(theta);
            varargout{1} = Xi;
        
        %%% Depth along the fault plane [km].
            % Depth = Xi*sin(theta);
    end
    
%%% Ice density [kg/m^3] as a function of temperature. From Fukusako (1990),
%%% Thermophysical Properties of Ice, Snow, and Sea Ice. There may be a more up to date
%%% reference for ice density. However, the temperature dependence does not have much
%%% affect relative to using a constant density.
    i_rho = TempC < -140;
    rho = 917*(1 - 1.17e-4*TempC);                              % eq. (4) in Fukusako
    rho(i_rho) = 930*(1 - 1.54e-5*TempC(i_rho));                % eq. (5) in Fukusako
    % rho = rho_ice*ones(numel(Depth),1);

%%% Lithostatic stress. Here we assume that sigma_3 = sigma_L.
    sigma_L = 1e-3*g*cumtrapz(Depth, rho);                      % [MPa]    

%%% Coefficient of friction.
    switch mu_flag
        case 'constant'
            mu = 0.55;                                      % Howell & Pappalardo, static
            % mu = 0.37;                                      % Howell & Pappalardo, kinetic

        case 'temperature'
        %%% Coefficient of friction from Persson (2015). 
        %%% 9.6.2024 - TempK_P is close to, but not exactly the same as TempK. So that 
        %%% needs to be addressed.
            [TempK_P, mu] = IceFriction(TempK);
    end

%%% Cohesion [MPa].
    C = 0;
    % C = 1;

%%% For a failure envelope, the angle of the frictional failure planes is determined from
%%% the friction coefficient.
    if strcmp(calc, 'failure')    
        theta = atan(mu);
    end

%%% Compute sigma_1 in the brittle regime using the friction law.
    sigma_1 = (2*C + (sin(2*theta) + mu*(cos(2*theta) + 1)).*sigma_L)...
        ./(sin(2*theta) + mu*(cos(2*theta) - 1));

%%% Compute the shear resistance or differential stress.
    switch calc
        
        case 'fault'        
        %%% Frictional shear resistance on the fault/failure plane in terms of principal 
        %%% stresses.
            Tau_F = 0.5*(sigma_1 - sigma_L)*sin(2*theta);
            
        case 'failure'        
        %%% For a failure enevelope, we define Tau_F as the differential stress:
        %%% sigma_1 - sigma_L. Since we are assuming that sigma_3 = sigma_L.
            Tau_F = sigma_1 - sigma_L;

    end

%%% Ductile differential stress from ice flow laws.
    Sigma_Diff = IceFlowLaws(TempK);
    
%%% Find the shear resistance.
    S = [Tau_F, Sigma_Diff];
    S = min(S,[],2);
    
    switch calc
        case 'fault'
        %%% Minimum force per unit distance)to promote incipient convergence is found by
        %%% integrating the shear resistance along the fault plane. This force is directed
        %%% along the fault plane, i.e. it is not the force in the sigma_1 (horizontal)
        %%% direction.
            F = 1e9*trapz(Xi, S);                                   % [N / m]

        case 'failure'
        %%% Force per unit distance in the sigma_1 direction (CHECK THIS) to induce failure.
            F = 1e9*trapz(Depth, S);                                % [N / m]
    end

end