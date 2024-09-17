function [Temp_P, mu_P] = IceFriction(TempK)
%%% This function carries out calculations to predict the temperature
%%% dependence of friction coeficient during steady sliding of ice on ice.

%%% Schuslon calculation.
    % [Temp, mu_S] = SchulsonIce;

%%% Persson calculation.
    Tc = 273.15;
    v = 1e-5;
    
    Temp_P = nan(numel(TempK), 1);
    mu_P = nan(numel(TempK), 1);
    for i = 1:numel(TempK)
        T_0 = TempK(i) - Tc;
        [T, mu_k] = PerssonIce(v, T_0);
        Temp_P(i, 1) = T;
        mu_P(i, 1) = mu_k;
    end
    Temp_P = Temp_P + Tc;

end

function [x, mu_k] = PerssonIce(v, T_0)

%%% Parameters from page 8.
l = 0.01;               % [m]
sigma_Y = 50e6;         % [Pa]
rho = 1e3;              % [kg / m^3]
kappa = 2;              % [W / (m K)]
c_p = 2e3;              % [J / (kg K)]
tau_m_0 = 50e6;         % [Pa]
beta = 0.33;
T_c = 273.15;           % [K]

%%% Properties of Barre granite.
% rho_g = 2640;
% kappa_g = 2.72;
% c_p_g = 904.4;

rho_g = rho;
kappa_g = kappa;
c_p_g = c_p;


%%% Compute temperature at the asperity surface;
%C = (l*v./(pi*rho*c_p*kappa)).^(1/2).*tau_m_0;

%alpha = (pi*rho*c_p*kappa*v./l).^(1/2);
alpha = (pi*v/l)^(1/2)*((rho*c_p*kappa)^(1/2) + (rho_g*c_p_g*kappa_g)^(1/2))/2;
C = v*tau_m_0./alpha;

fun = @(T)(T - (T_c + T_0) - C.*(1 - T./T_c).^beta);
%%% Need to automatically set T_init lower bound.
T_int = [-230 0] + T_c;
x = fzero(fun, T_int) - T_c;

%%% Compute shear stress on asperity;
tau_m = alpha.*(x - T_0)./v;

%%% Compute friction coefficient.
%sigma_Y = 1e6*(-5.08*x + 15.19);
mu_k = tau_m./sigma_Y;

end

function [Temp, mu_kw] = SchulsonIce

R = 8.3144626181;             % [J/(mol K)] gas constant
Tc = 273.15;
% T = Tc + linspace(-20, -0.5, 100)';  % [K] temperature range
T = linspace(100, Tc - 0.5, 100)';
Temp = T - Tc;

%%%-------------------------------------------------------------------------------%%%
%%% Parameters for Schulson 2015 calculation of ice-ice friction coefficient, equation (38).
%%% Warm ice, bottom of page 20
  Q1 = 145e3; % [J/mol] indentation activation energy
  Q2 = 120e3; % [J/mol] creep rate shear stress activation energy
  n1 = 4.2;
  n2 = 3;
  
%%% "creep rate" shear stress used by Schulson 2015.
  v_s = 1e-5;          % [m/s] sliding velocity
  h = 1e-4/2;            % [m] near-surface inelastic zone
  a = 15e-6;           % [m] asperity radius
  C2 = 0.0094;         % C2 = 1/(B*C), estimated from Schulson after equation (41)
  strainrate = v_s/h;
  time = v_s/(2*a);
  
  mu_k = C2*(strainrate.^(1/n2).*time.^(-1/n1)).*exp(Q2./(R*n2*T) - Q1./(R*n1*T));

%%% Warm ice from Schulson
  L_v = 320e6; % [J/m^3]
  kappa = 2.1; % [W/(m K)]
  rho = 917;   % [kg/m^3]
  c_p = 1.9e3; % [J/(kg K)]
  
%%% Properties of Barre granite.
  % rho_g = 2640;
  % kappa_g = 2.72;
  % c_p_g = 904.4;
  rho_g = rho;
  kappa_g = kappa;
  c_p_g = c_p;
 
  f = 100;
  
  delta = 1e-6;
  %dT = 10; % this is the value Schulson used in his example. It's temp below melting?
  dT = abs(T - Tc);
  
  %%% Equation (46)
  %t_c = (L_v*delta./(2*dT)).^2.*(1./(kappa*rho*c_p));
  t_c = (L_v*delta./dT).^2.*((kappa*rho*c_p)^(1/2) + (kappa_g*rho_g*c_p_g)^(1/2))^(-2);
  
  v_t = 2*a./(f*t_c);
  %v_t = 8*kappa*rho*c_p*a*dT.^2/(f*L_v^2*delta^2);
  
  gamma = 0.4;
  eta = gamma*log(v_s./v_t);
  mu_kw = (1 - eta).*mu_k;
end