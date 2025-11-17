function [S_Compression, Tau_Compression, S_Tension, Tau_Tension, Sigma_Diff]...
    = IceFailureEnvelope(Depth, TempK, sigma_L, mu_flag, varargin)

%%% This function computes a failure envelope by assuming:
%%% 1. frictional failure takes place on optimally oriented planes.
%%% 2. Strain rates for ductile failure stress are determined kinematically from the
%%%    slab centerline curvature, following Buffett (2006).

%%% 9.6.2024 - TempK_P is close to, but not exactly the same as TempK. So that needs 
%%% to be addressed. Also, the Persson friction values need to be tweaked a bit to 
%%% better fit the experimental data.

%%% 10.29.2024 - Need to differentiate between tension and compression.

%%% 10.29.2024 - Need to account for thin plate (plane strain?) stress state 
%%%              assumption: sigma_zz vanishes.

%%%-------------------------------------------------------------------------------%%%
%%% INPUT
%%%-------------------------------------------------------------------------------%%%
% Depth   - A vector of depth values (in meters) through the desired shell thickness.

% TempK   - A vector of temperature values (in Kelvin) for each value of depth.

% mu_flag - Determines how to handle the friction coefficient. 
%           Set mu_flag = 'temperature' to use temperature dependent values from 
%           Persson (2015).
%           Set mu_flag = 'constant' to use a constant value.


%%%-------------------------------------------------------------------------------%%%
%%% OUTPUT
%%%-------------------------------------------------------------------------------%%%
% S       - The failure envelope.

% F       - The integrated force with depth.
%%%-------------------------------------------------------------------------------%%%

    if isempty(varargin) == false
        dedt = varargin{1};                                 % strain rate [1/s]
    else
        dedt = 1e-15;
    end

%%% Ice density [kg/m^3] as a function of temperature. From Fukusako (1990),
%%% Thermophysical Properties of Ice, Snow, and Sea Ice. There may be a more up to date
%%% reference for ice density. However, the temperature dependence does not have much
%%% affect relative to using a constant density.
    % i_rho = TempC < -140;
    % rho = 917*(1 - 1.17e-4*TempC);                              % eq. (4) in Fukusako
    % rho(i_rho) = 930*(1 - 1.54e-5*TempC(i_rho));                % eq. (5) in Fukusako
    % rho = rho_ice*ones(numel(Depth),1);

%%% Lithostatic stress. Here we assume that sigma_3 = sigma_L.
    % sigma_L = 1e-6*g*cumtrapz(Depth, rho);                      % [MPa]
    % sigma_L = 1e-6*g*(cumtrapz(Depth, rho) + mean(rho)*Depth(1));                      % [MPa]

%%% Coefficient of friction.
    switch mu_flag
        case 'constant'
            % mu = 0.55;                                % Howell & Pappalardo, static
            mu = 0.37;                                % Howell & Pappalardo, kinetic

        case 'temperature'
        %%% Coefficient of friction from Persson (2015).        
            [TempK_P, mu] = IceFriction(TempK);
    end

%%% Cohesion [MPa].
    C = 0;
    % C = 1;

%%% For a failure envelope, the angle of the frictional failure planes is determined 
%%% from the friction coefficient.
    phi = atan(mu);

%%%-------------------------------------------------------------------------------%%%
%%% Compression.
%%%-------------------------------------------------------------------------------%%%
%%% For compression: sigma_3 = sigma_L.
    sigma_3 = sigma_L;

%%% Compute sigma_1 for compression in the brittle regime using the friction law. 
%%% Brice & Kohlstedt equation.
    sigma_1 = (2*C + (cos(phi) + mu.*(sin(phi) + 1)).*sigma_3)./...
        (cos(phi) + mu.*(sin(phi) - 1));
    % sigma_1 = (2*C + (sin(2*phi) + mu.*(cos(2*phi) + 1)).*sigma_L)...
    %     ./(sin(2*phi) + mu.*(cos(2*phi) - 1));

%%% Differential stress at failure for compression.
    Tau_Compression = sigma_1 - sigma_3;

%%% Ductile differential stress from ice flow laws. Remove complex values that result 
%%% from negative strain rates (i.e. tension).
    Sigma_Diff = IceFlowLaws(TempK, dedt);
    Sigma_Diff(imag(Sigma_Diff) ~= 0) = nan;
    
%%% Find the failure envelope for compression.
    S_Compression = [Tau_Compression, Sigma_Diff];
    S_Compression = min(S_Compression,[],2);

%%%-------------------------------------------------------------------------------%%%
%%% Tension.
%%%-------------------------------------------------------------------------------%%%
%%% For tension: sigma_1 = sigma_L and sigma_3 = sigma_H.
    sigma_1 = sigma_L;

%%% Sigma_3 for tension below 200 MPa normal stress.
    sigma_3 = ((cos(phi) + mu.*(sin(phi) - 1)).*sigma_1 - 2*C)./...
        (cos(phi) + mu.*(sin(phi) + 1));
    
%%% Differential stresses at failure for tension.
    Tau_Tension = sigma_1 - sigma_3;

%%% Find the failure envelope for tension.
    S_Tension = [Tau_Tension, Sigma_Diff];
    S_Tension = -min(S_Tension,[],2);
       
%%% Convert to Pa.
    S_Compression = 1e6*S_Compression;
    Tau_Compression = 1e6*Tau_Compression;
    S_Tension = 1e6*S_Tension;
    Tau_Tension = 1e6*Tau_Tension;
    Sigma_Diff = 1e6*Sigma_Diff;

end