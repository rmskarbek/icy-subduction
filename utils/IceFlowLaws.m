function Sigma_Diff = IceFlowLaws(TempK, varargin)

%%% This function evaluates the ductile differential stress using flow laws for ice.
%%% This function is called by: IceFailureEnvelope.m

%%%-------------------------------------------------------------------------------%%%
%%% INPUT
%%%-------------------------------------------------------------------------------%%%
% TempK    - A vector of temperature values (in Kelvin) for each value of depth.

% varargin - Use to input strain rate as a function of temperature. Otherwise strain 
%            rate is set to 1e-15 [1/s]

%%%-------------------------------------------------------------------------------%%%
%%% OUTPUT
%%%-------------------------------------------------------------------------------%%%
% Sigma_Diff - A vector of differential stress values corresponding to the input
%              temperature values.

%%%-------------------------------------------------------------------------------%%%
    R = 8.314;                                         % gas constant [J / (mol K)]

%%% Set the strain rate if it is not input.
    if isempty(varargin)
        StrainRate = 1e-15;                            % strain rate [1/s]
    else
        StrainRate = abs(varargin{1});                 % strain rate [1/s]
    end
    
    
%%% Indices for high temperature behavior.
    T_tr = 259;                                        % low/high transition temp [K]
    i_low = TempK <= T_tr;

%%%-------------------------------------------------------------------------------%%%
%%% Dislocation creep from Behn et al., (2021). Low temperatures are less than or 
%%% equal to 259 K.
%%%-------------------------------------------------------------------------------%%%
    A_disl_lowT = 4e4;                                      % [1 / (s MPa^n m^m)]
    A_disl_highT = 6e28;                                    % [1 / (s MPa^n m^m)]
    m_disl = 0;                                             % grain size exponent
    n_disl = 4;                                             % stress exponent
    Q_disl_lowT = 60e3;                                     % [J / mol]
    Q_disl_highT = 180e3;                                   % [J / mol]

%%% Differential stress for low temp dislocation creep [MPa].
    sigma_disl_lowT = ((StrainRate/A_disl_lowT)...
        .*exp(Q_disl_lowT./(R*TempK))).^(1/n_disl);

%%% Composite differential stress for dislocation creep. First compute the high 
%%% temperature stress, then replace appropriate elements with the low temperature 
%%% stress.
    sigma_disl = ((StrainRate/A_disl_highT)...
        .*exp(Q_disl_highT./(R*TempK))).^(1/n_disl);
    sigma_disl(i_low) = sigma_disl_lowT(i_low);

%%% Conversion factor from Dombard & McKinnon (2001).
    A_prime_disl = (3^(1/2)/2)^(n_disl + 1);
    sigma_disl = sigma_disl/A_prime_disl;


%%%-------------------------------------------------------------------------------%%%
%%% Grain boundary sliding from Behn et al., (2021). Low temperatures are less than 
%%% or equal to 259 K.
%%%-------------------------------------------------------------------------------%%%
    d = 10e-6;                                              % grain size [m]
    A_gbs_lowT = 3.9e-3;                                    % [1 / (s MPa^n m^m)]
    A_gbs_highT = 3e26;                                     % [1 / (s MPa^n m^m)]
    m_gbs = 1.4;                                            % grain size exponent
    n_gbs = 1.8;                                            % stress exponent
    Q_gbs_lowT = 49e3;                                      % [J / mol]
    Q_gbs_highT = 192e3;                                    % [J / mol]

%%% Differential stress for low temp grain boundary sliding [MPa].
    sigma_gbs_lowT = ((StrainRate*d^m_gbs/A_gbs_lowT)...
        .*exp(Q_gbs_lowT./(R*TempK))).^(1/n_gbs);

%%% Composite differential stress for grain boundary sliding. First compute the high 
%%% temperature stress, then replace appropriate elements with the low temperature 
%%% stress.
    sigma_gbs = ((StrainRate*d^m_gbs/A_gbs_highT)...
        .*exp(Q_gbs_highT./(R*TempK))).^(1/n_gbs);
    sigma_gbs(i_low) = sigma_gbs_lowT(i_low);

%%% Conversion factor from Dombard & McKinnon (2001).
    A_prime_gbs = (3^(1/2)/2)^(n_gbs + 1);
    sigma_gbs = sigma_gbs/A_prime_gbs;


%%%-------------------------------------------------------------------------------%%%
%%% Basal slip-accommodated GBS from Goldsy & Kohlstedt (2001), Table 5. From 
%%% Christine and Behn et al. (2021) this mechanism is not expected to be active at 
%%% grain sizes expected in terrestrial settings (d > 1 mm), but could be active at 
%%% grain sizes on icy moons where it is very cold and grain growth is suppressed.
%%%-------------------------------------------------------------------------------%%%
    A_bs = 5.5e7;                                           % [1 / (s MPa^n m^m)]
    m_bs = 0;                                               % grain size exponent
    n_bs = 2.4;                                             % stress exponent
    Q_bs = 60e3;                                            % [J / mol]

%%% Differential stress for basal slip-accommodated GBS [MPa].
    sigma_bs = ((StrainRate/A_bs).*exp(Q_bs./(R*TempK))).^(1/n_bs);

%%% Conversion factor from Dombard & McKinnon (2001).
    A_prime_bs = (3^(1/2)/2)^(n_bs + 1);
    sigma_bs = sigma_bs/A_prime_bs;


%%%-------------------------------------------------------------------------------%%%
%%%-------------------------------------------------------------------------------%%%
%%% Combine the mechanisms.
    Sigma_Diff = sigma_disl + sigma_gbs + sigma_bs;

end