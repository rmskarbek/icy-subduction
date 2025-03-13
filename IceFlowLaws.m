function Sigma_Diff = IceFlowLaws(TempK, varargin)

% varargin - Use to input strain rate as a function of temperature. Otherwise strain rate
%            is set to 1e-15 [1/s]

%%% 9.30.2024 - Add grain size as an input? Currently it is set on Line 39.

%%% 10.22.2024 - ductile mechanisms act in parallel.

    R = 8.314;                                              % gas constant [J / (mol K)]

%%% Set the strain rate if it is not input.
    if isempty(varargin)
        e = 1e-15;                                          % strain rate [1/s]
    else
        % e = varargin{1};                                    % strain rate [1/s]
        e = abs(varargin{1});                                    % strain rate [1/s]
    end
    
    
%%% Indices for high temperature behavior.
    T_tr = 259;                                             % low/high transition temp [K]
    i_low = TempK <= T_tr;

%----------------------------------------------------------------------------------------%
%----------------------------------------------------------------------------------------%
%%% Dislocation creep from Behn et al., (2021). Low temperatures are less than or equal to
%%% 259 K.
    A_disl_lowT = 4e4;                                      % [1 / (s MPa^n m^m)]
    A_disl_highT = 6e28;                                    % [1 / (s MPa^n m^m)]
    m_disl = 0;                                             % grain size exponent
    n_disl = 4;                                             % stress exponent
    Q_disl_lowT = 60e3;                                     % [J / mol]
    Q_disl_highT = 180e3;                                   % [J / mol]

%%% Differential stress for low temp dislocation creep [MPa].
    sigma_disl_lowT = ((e/A_disl_lowT).*exp(Q_disl_lowT./(R*TempK))).^(1/n_disl);

%%% Composit differential stress for dislocation creep. First compute the high temperature
%%% stress, then replace appropriate elements with the low temperature stress.
    sigma_disl = ((e/A_disl_highT).*exp(Q_disl_highT./(R*TempK))).^(1/n_disl);
    sigma_disl(i_low) = sigma_disl_lowT(i_low);

%%% Conversion factor from Dombard & McKinnon (2001).
    A_prime_disl = (3^(1/2)/2)^(n_disl + 1);
    sigma_disl = sigma_disl/A_prime_disl;

%----------------------------------------------------------------------------------------%
%----------------------------------------------------------------------------------------%
%%% Grain boundary sliding from Behn et al., (2021). Low temperatures are less than or 
%%% equal to 259 K.
    d = 10e-6;                                              % grain size [m]
    A_gbs_lowT = 3.9e-3;                                    % [1 / (s MPa^n m^m)]
    A_gbs_highT = 3e26;                                     % [1 / (s MPa^n m^m)]
    m_gbs = 1.4;                                            % grain size exponent
    n_gbs = 1.8;                                            % stress exponent
    Q_gbs_lowT = 49e3;                                      % [J / mol]
    Q_gbs_highT = 192e3;                                    % [J / mol]

%%% Differential stress for low temp grain boundary sliding [MPa].
    sigma_gbs_lowT = ((e*d^m_gbs/A_gbs_lowT).*exp(Q_gbs_lowT./(R*TempK))).^(1/n_gbs);

%%% Composit differential stress for grain boundary sliding. First compute the high 
%%% temperature stress, then replace appropriate elements with the low temperature stress.
    sigma_gbs = ((e*d^m_gbs/A_gbs_highT).*exp(Q_gbs_highT./(R*TempK))).^(1/n_gbs);
    sigma_gbs(i_low) = sigma_gbs_lowT(i_low);

%%% Conversion factor from Dombard & McKinnon (2001).
    A_prime_gbs = (3^(1/2)/2)^(n_gbs + 1);
    sigma_gbs = sigma_gbs/A_prime_gbs;

%----------------------------------------------------------------------------------------%
%----------------------------------------------------------------------------------------%
%%% Basal slip-accommodated GBS from Goldsy & Kohlstedt (2001), Table 5. From Christine
%%% and Behn et al. (2021) this mechanism is not expected to be active at grain sizes
%%% expected in terrestrial settings (d > 1 mm), but could be active at grain sizes on
%%% icy moons where it is very cold and grain growth is suppressed.
    A_bs = 5.5e7;                                           % [1 / (s MPa^n m^m)]
    m_bs = 0;                                               % grain size exponent
    n_bs = 2.4;                                             % stress exponent
    Q_bs = 60e3;                                            % [J / mol]

%%% Differential stress for basal slip-accommodated GBS [MPa].
    sigma_bs = ((e/A_bs).*exp(Q_bs./(R*TempK))).^(1/n_bs);

%%% Conversion factor from Dombard & McKinnon (2001).
    A_prime_bs = (3^(1/2)/2)^(n_bs + 1);
    sigma_bs = sigma_bs/A_prime_bs;

%----------------------------------------------------------------------------------------%
%----------------------------------------------------------------------------------------%
%%% Combine the mechanisms.
    % Sigma_Diff = [sigma_disl, sigma_gbs, sigma_bs];
    Sigma_Diff = sigma_disl + sigma_gbs + sigma_bs;


%----------------------------------------------------------------------------------------%
%----------------------------------------------------------------------------------------%
%%% Ice flow laws, paramters for Regimes A, B, and C are from from Dombard & McKinnon 
%%% (2001), Table 1. These regimes represent grainsize insensitive disclocation creep 
%%% mechanisms. See Durham & Stern (2001) for references. See text in D & M (2006) in the
%%% paragraph after eq (1b). Symbols here are after Dombard & McKinnon (2001).

%%% Regime A, 240 - 258 K.
    % Q_A = 91e3;                                             % [J / mol]
    % n_A = 4;                                                % stress exponent
    % A_A = 10^11.8;                                          % [1 / (s MPa^n)]
    % A_A = A_A*(3^(1/2)/2)^(n_A + 1);
    % Tau_A = ((e/A_A)*exp(Q_A./(R*TempK))).^(1/n_A);         % [MPa]
    
%%% Regime B.
    % Q_B = 61e3;                                             % [J / mol]
    % n_B = 4;                                                % stress exponent
    % A_B = 10^5.1;                                           % [1 / (s MPa^n)]
    % A_B = A_B*(3^(1/2)/2)^(n_B + 1);
    % Tau_B = ((e/A_B)*exp(Q_B./(R*TempK))).^(1/n_B);         % [MPa]
    
%%% Regime C.   
    % Q_C = 39e3;                                             % [J / mol]
    % n_C = 6;                                                % stress exponent
    % A_C = 10^(-3.8);                                        % [1 / (s MPa^n)]
    % A_C = A_C*(3^(1/2)/2)^(n_C + 1);
    % Tau_C = ((e/A_C)*exp(Q_C./(R*TempK))).^(1/n_C);         % [MPa]
    
%%% Grainsize sensitive grain boundary sliding from Goldsby & Kohlstedt (2001). Here,
%%% using parameter values listed in Table in in Dombard & McKinnon (2001).
    % d = 1e-3;                                               % grainsize [m]
    % Q_G = 49e3;                                             % [J / mol]
    % m_G = 1.4;                                              % grainsize exponent
    % n_G = 1.8;                                              % stress exponent
    % A_G = 3.9e-3;                                           % [1 / (s m^m MPa^n)]
    % A_G = A_G*(3^(1/2)/2)^(n_G + 1);
    % Tau_G = ((e*d^m_G/A_G)*exp(Q_G./(R*TempK))).^(1/n_G);         % [MPa]

    % Sigma_Diff = [Tau_A, Tau_B, Tau_C, Tau_G];
end
