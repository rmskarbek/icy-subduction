function rho_ice = IceDensity(TempK, method)

%%% This function computes the density of ice for a given temperature, using 
%%% equations (4) and (5) in:
%%%     Fukusako 1990 - Thermophysical properties of ice, snow, and sea ice

%%%-------------------------------------------------------------------------------%%%
%%% Notes.
%%%-------------------------------------------------------------------------------%%%
%%% 1. Both Johnson2017 and Howell2019 use a temperature-dependent coefficient of 
%%%    volumetric thermal expansion to calculate the density of solid ice. They use 
%%%    an equation referenced to Kirk & Stevenson (1987), who in turn cite 
%%%    Hobbs (1974).

%%% 2. Here I use an equation for the ice density from Fukusako (1990) that is based 
%%%    on the data from Hobbs and some other papers.

%%% 3. The thermal expansion method results in a maximum ice density at ~150 K that 
%%%    does not agree with the density data presented in Figure 4 of Fukusako. This 
%%%    density behavior is embedded in the codes used by Johnson2017 and Howell2019.


%%%-------------------------------------------------------------------------------%%%
%%% INPUT
%%%-------------------------------------------------------------------------------%%%
% TempK  - A vector of temperature values (in Kelvin) for each value of depth.

% method - A string that determines how to compute the density.
%          method = 'constant' - uses a constant density.

%          method = 'density'  - uses equation for the ice density as a function of
%          temperature from Fukusako (1990).

%          method = 'thermal' - uses a coefficient of thermal expansion to compute the
%          density as a function of temperature.

%%%-------------------------------------------------------------------------------%%%
%%% OUTPUT
%%%-------------------------------------------------------------------------------%%%
% rho_ice - A vector of ice density values corresponding to the input temperature
%           values.

%%%-------------------------------------------------------------------------------%%%
Tc = 273.15;                                                % water melt temp 1 atm [K]

switch method
    case 'constant'
        rho_ice = 917;
        
    case 'density'
        TempC = TempK - Tc;                                 % [C]
        T_transition = -140;
        i_rho = TempC < T_transition;
        
%%% Eq. (4) in Fukusako
        rho_ice = 917*(1 - 1.17e-4*TempC);                  % [kg/m^3]

%%% Eq. (5) in Fukusako
        rho_ice(i_rho) = 930*(1 - 1.54e-5*TempC(i_rho));    % [kg/m^3]

    case 'thermal'
        % alpha = 1e-6*(0.67*TempK - 24.86);                  % [1/K] eq. (6) in Fukusako

%%% Coefficient of thermal expansopn via Hobbs? See Johnson2017, bottom p. 2770
        alpha = (1.56e-4/250)*TempK;                        % [1/K]
        rho_ice = 917*(1 + alpha.*(Tc - TempK));
end