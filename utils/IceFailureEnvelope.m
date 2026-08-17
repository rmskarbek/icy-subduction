function [S_Compression, Tau_Compression, S_Tension, Tau_Tension, Sigma_Diff]...
    = IceFailureEnvelope(Depth, TempK, sigma_L, mu_flag, varargin)

%%% This function computes a failure envelope as a function of temperature and
%%% pressure by evaluated models for both brittle and ductile failure. The failure
%%% envelope stress is defined by whichever mechanism yields a smaller failure stress.
%%% This function is called by BendingForce.m

%%%-------------------------------------------------------------------------------%%%
%%% INPUT
%%%-------------------------------------------------------------------------------%%%
% Depth    - A vector of depth values (in meters) through the desired shell thickness.

% TempK    - A vector of temperature values (in Kelvin) for each value of depth.

% Sigma_L  - A vector of lithostatic stress values for each value of depth.

% mu_flag  - A string that determines how to handle the friction coefficient. 
%            Set mu_flag = 'temperature' to use temperature dependent values from 
%            Persson (2015).
%            Set mu_flag = 'constant' to use a constant value.

% varargin - Use to input strain rate as a function of depth. Otherwise strain 
%            rate is set to 1e-15 [1/s]

%%%-------------------------------------------------------------------------------%%%
%%% OUTPUT
%%%-------------------------------------------------------------------------------%%%
% S_Compression   - The failure envelope for compression.

% S_Tension       - The failure envelope for tension.

% Tau_Compression - Frictional failure stress for compression.

% Tau_Tension     - Frictional failure stress for tension.

% Sigma_Diff      - Ductile failure stress from ice flow laws, computed with 
%                   IceFlowLaws.m

%%%-------------------------------------------------------------------------------%%%

    if isempty(varargin) == false
        dedt = varargin{1};                                 % strain rate [1/s]
    else
        dedt = 1e-15;
    end

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

%%% For a failure envelope, the angle of the frictional failure planes is determined 
%%% from the friction coefficient.
    phi = atan(mu);

%%%-------------------------------------------------------------------------------%%%
%%% Compression.
%%%-------------------------------------------------------------------------------%%%
%%% For compression: sigma_3 = sigma_L.
    sigma_3 = sigma_L;

%%% Compute sigma_1 for compression in the brittle regime using the friction law. 
%%% Brace & Kohlstedt equation.
    sigma_1 = (2*C + (cos(phi) + mu.*(sin(phi) + 1)).*sigma_3)./...
        (cos(phi) + mu.*(sin(phi) - 1));

%%% Differential stress at failure for compression.
    Tau_Compression = sigma_1 - sigma_3;

%%% Ductile differential stress from ice flow laws. Remove complex values that result 
%%% from negative strain rates (i.e. tension).
    Sigma_Diff = IceFlowLaws(TempK, dedt);
    Sigma_Diff(imag(Sigma_Diff) ~= 0) = nan;
    
%%% Find the failure envelope for compression.
    % Rheol_Compression = ones(size(TempK));
    S_Compression = [Tau_Compression, Sigma_Diff];
    % Rheol_Compression(Sigma_Diff < Tau_Compression) = 2;
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
    % Rheol_Tension = ones(size(TempK));
    S_Tension = [Tau_Tension, Sigma_Diff];
    % Rheol_Tension(Sigma_Diff < Tau_Tension) = 2;
    S_Tension = -min(S_Tension,[],2);
       
%%% Convert to Pa.
    S_Compression = 1e6*S_Compression;
    Tau_Compression = 1e6*Tau_Compression;
    S_Tension = 1e6*S_Tension;
    Tau_Tension = 1e6*Tau_Tension;
    Sigma_Diff = 1e6*Sigma_Diff;

end