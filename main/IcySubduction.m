function Out = IcySubduction(p)

%%% This code solves the dimensional equations for temperature and porosity evolution 
%%% in a subducting 1-D column of ice. The code is essentially the same model used by
%%% Johnson et al. (2017) and Howell & Papallardo (2019).

%%%-------------------------------------------------------------------------------%%%
%%% INPUT
%%%-------------------------------------------------------------------------------%%%
% p   - A matlab structure containting parameter values for the simulation. Generate 
%       with: p = IcySubduction_Parameters.m


%%%-------------------------------------------------------------------------------%%%
%%% OUTPUT. All of these variables are stored in a structure called 'Out'.
%%%-------------------------------------------------------------------------------%%%
% Out - A structure that contains the simulation output, the slab geometry, and also 
%       the same parameters structure (p) that is input.

% p              - the same parameters structure that is input.

% Density        - [kg/m^3] Density in the column at each time step.

% Porosity       - [-] Porosity in the column at each time step.

% Temperature    - [K] Temperature in the column at each time step.

% TimeYrs        - [years] The time at each step that the solver routine used.

% VerticalStress - [Pa] Vertical stress (i.e. overburden stress) in the column at
%                  each time step.

% Viscosity      - [Pa*s] Viscosity in the column at each time step.

%%%-------------------------------------------------------------------------------%%%
%%% The following variables are stored in a structure called 'Out.Slab', that 
%%% contains geometric information of the subducting slab.

% ArcLength      - [m] Length along the center line of the slab for each time step.

% Depth          - [m] Depth in the column for each time step, measured from the
%                  upper surface of the ice shell.

% CenterLine     - [m] Cartesian coordinates of the centerline for each time step,
%                  stored as a complex number: z = x + iy.

% Top            - [m] Cartesian coordinates of the upper surface of the slab for 
%                  each time step, stored as a complex number: z = x + iy.

% Bottom         - [m] Cartesian coordinates of the upper surface of the slab for 
%                  each time step, stored as a complex number: z = x + iy.

% Dip            - [radians] Dip angle of the center line for each time step.

% Curvature      - [1 / m] Curvature of the center line for each time step.

% dKds           - [1 / m^2] Curvature gradient of the center line for each time
%                  step, computed along the arc length coordinate (s).


%%%-------------------------------------------------------------------------------%%%
%%% Constants.
%%%-------------------------------------------------------------------------------%%%
spy = 365.25*24*3600;                   % seconds per year
g = p.Constants.Gravity;                % Europa gravity [m/s^2]
R = p.Constants.GasConst;               % gas constant [J/(mol K)]


%%%-------------------------------------------------------------------------------%%%
%%% Geometry
%%%-------------------------------------------------------------------------------%%%
H = p.Geometry.SlabThick;               % thickness of conductive slab [m]
H_shell = p.Geometry.ShellThick;        % total thickness of ice shell [m]
GeoFlag = p.Geometry.GeoFlag;           % a flag that determines the slab geometry.
R_min = p.Geometry.CurveRadius;         % minimum radius of plate curvature [m]
v_plate = p.Geometry.PlateRate;         % plate convergence rate [m/s]
% s_end = p.Geometry.ArcLength;           % total arc length of center line [m]
type = 'location';                      % flag for Buffett2006.m

if strcmp(p.Geometry.GeoFlag, 'CircularArc')
    DipAngle = p.Geometry.DipAngle;     % final slab dip angle [degrees]
end

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
f_slab = p.Porosity.SaltSlab;           % salt content in slab
f_conduct = p.Porosity.SaltConduct;     % salt content in non-subducting conductive layer
f_convect = p.Porosity.SaltConvect;     % salt content in convective layer

method = p.Porosity.DensityMethod;      % method for computing ice density
rho_ice_b = IceDensity(T_b, method);    % density of convecting ice.
rho_salt = p.Porosity.SaltDensity;      % salt density [kg/m^3]
rho_ice_conv = (1 - f_convect)*rho_ice_b + f_convect*rho_salt;

%%%-------------------------------------------------------------------------------%%%
%%% Numerical
%%%-------------------------------------------------------------------------------%%%
N = p.Numerical.GridPoints;             % Number of grid points
dz = p.Numerical.GridSpacing;           % Numerical grid spacing [m]

%%% Flag for how to handle variable thermal conductivity in the finite difference 
%%% scheme.
coeff = p.Numerical.VarCoeff;

RunTime = p.Numerical.RunTime;          % run time of the simulation [s]
time = linspace(0, RunTime, 1e3);       % vector of time values for solution output.


%%%-------------------------------------------------------------------------------%%%
%%% Initial conditions.
%%%-------------------------------------------------------------------------------%%%
%%% Initial conditions.
Vars_init = p.Numerical.InitialCond;

%%% Depth referenced to the top of the subducting column.
Depth_C = linspace(0, H, N)';

%%% Initial porosity, ice density, and bulk density in the slab.
phi_init = p.Numerical.InitialCond(N+1:2*N);
rho_ice_init = IceDensity(Vars_init(1:N,1), method);

%%% Lithostatic stress in non-subducting conductive portion of the ice shell. This 
%%% assumes the same porosity profile as in the slab.
rho_cond = ((1 - f_conduct)*rho_ice_init + f_conduct*rho_salt).*(1 - phi_init);
stress_cond = cumtrapz(Depth_C, g*rho_cond);

    
%%%-------------------------------------------------------------------------------%%%
%%% Conduct the simulation and deal with the output.
%%%-------------------------------------------------------------------------------%%%
%%% Generate the sparsity pattern for amazing speed up!
S = SparsityPattern(N);
options = odeset('JPattern', S, 'RelTol', 1e-10,'AbsTol', [1e-12*ones(N,1);...
    1e-12*ones(N,1)]);
% options = odeset('JPattern', S);
% options = odeset('JPattern', S, 'Vectorized', 'on');

%%% Solve the system.
[Time, Vars_f] = ode15s(@PDE, time, Vars_init, options);
Temp = Vars_f(:, 1:N)';
Phi = Vars_f(:, N+1:2*N)';
ArcLength = Time*v_plate;                                                       % [m]

%%% Compute final ice and bulk densities in the slab.
Rho_ice = IceDensity(Temp, method);
Rho = ((1 - f_slab)*Rho_ice + f_slab*rho_salt).*(1 - Phi);

%%% Evalute the geometry for the output time.
switch GeoFlag
    case 'Buffett'
        [Z_Top, Z_Center, Z_Bottom, Theta, K, dKds] = Buffett2006(type, H, H_shell,...
            R_min, ArcLength);

    case 'CircularArc'
        K = 1/R_min;
        dKds = 0;

        [Z_Top, Z_Center, Z_Bottom, Theta] = CircularArc(type, H, H_shell, R_min,...
            DipAngle, ArcLength);
end

%%% Get the distance and depth of the slab column at each output time and compute the 
%%% lithostatic stress (i.e. sigma_yy).
Distance = nan(N, numel(Time));
Depth = nan(N, numel(Time));
Sigma_yy = nan(N, numel(Time));

%%% Lithostatic stress referenced to the top of the subducting column.
    S_column = cumtrapz(Depth_C, g*Rho);

for i = 1:numel(Time)
    Distance(:, i) = linspace(real(Z_Top(i,1)), real(Z_Bottom(i,1)), N)';
    Depth(:, i) = linspace(imag(Z_Top(i,1)), imag(Z_Bottom(i,1)), N)';

%%% Total lithostatic stress in the column.
    if Depth(1,i) <= H
        Sigma_yy(:,i) = S_column(:,i) + interp1(Depth_C, stress_cond, Depth(1,i));
    else

%%% Vertical stress due to thickness of convecting ice above the top of the column.
        S_conv = g*rho_ice_conv.*(Depth(1,i) - H);

%%% Total lithostatic stress.
        Sigma_yy(:,i) = S_column(:,i) + stress_cond(end) + S_conv;
    end
end

%%% Convert time to years.
Time = Time/spy;

%%% Evaluate the viscosity in the slab.
Eta = eta_b*exp((Q/R)*(1./Temp - 1/T_b));
Eta(Eta > eta_max) = eta_max;

%%% Create output structure.
Out = struct('TimeYrs', Time, 'Temperature', Temp, 'Porosity', Phi, 'Density', Rho,...
    'VerticalStress', Sigma_yy, 'Viscosity', Eta, 'Slab', [], 'p', p);
Out.Slab.ArcLength = ArcLength;
Out.Slab.Distance = Distance;
Out.Slab.Depth = Depth;
Out.Slab.CenterLine = Z_Center;
Out.Slab.Top = Z_Top;
Out.Slab.Bottom = Z_Bottom;
Out.Slab.Dip = Theta;
Out.Slab.Curvature = K;
Out.Slab.dKds = dKds;


%%%-------------------------------------------------------------------------------%%%
%%% Numerical evolution equations.
%%%-------------------------------------------------------------------------------%%%
function dVarsdt = PDE(time, Vars, ~)
    T = Vars(1:N,:);
    phi = Vars(N+1:2*N,:);

%%% Get the current distance and burial depth of the numerical grid points.
    s = v_plate*time;
    switch GeoFlag
        case 'Buffett'
            [z_top, ~, z_bottom, theta] = Buffett2006(type, H, H_shell, R_min, s);

        case 'CircularArc'
            [z_top, ~, z_bottom] = CircularArc(type, H, H_shell, R_min, DipAngle, s);
    end
    % distance = linspace(real(z_top), real(z_bottom), N)';
    depth = linspace(imag(z_top), imag(z_bottom), N)';

%%% Compute ice density and bulk density in the slab for current
%%% temperature and porosity.
    rho_ice = IceDensity(T, method);
    rho = ((1 - f_slab)*rho_ice + f_slab*rho_salt).*(1 - phi);

%%%-------------------------------------------------------------------------------%%%
%%% Lithostatic (vertical) stress in the column.

%%% Compute the height of slab material directly above each location in the column.
%%% Buffett2006_Top.m will do this calculation exactly, but slows down the code a
%%% lot, ~200 times longer to excute.
    % y_top = Buffett2006_Top(s_end, H, H_shell, R_min, distance);
    % height_slab = depth - y_top;
    
%%% Instead, the height of the slab material can be approximated using the local slab
%%% dip angle. This does not slow down the code at all.
    height_slab = Depth_C./cos(theta);
    y_top = depth - height_slab;

%%% Lithoststic stress due to the height of slab above each location in the column.
%%% Here we have to make an assumption about the density of the overlaying slab 
%%% material. The most accurate thing would be to use the output of previous
%%% timesteps to interpolate the slab density along vertical profiles above each
%%% column location. This would be difficult and would certainly slow the code down
%%% by orders of mangitude. Instead, we just assume that the density of overlying
%%% material is the same as the column material.
    s_column = cumtrapz(height_slab, g*rho);

%%% Total lithostatic stress in the column.
    if depth(1,1) <= H
        s_lith = s_column + interp1(Depth_C, stress_cond, depth(1,1));
    else

%%% Vertical stress due to height of convecting ice above each location in the column.
        height_conv = y_top - H;
        s_conv = g*rho_ice_conv.*height_conv;

%%% Total lithostatic stress.
        s_lith = s_column + stress_cond(end) + s_conv;
    end

%%%-------------------------------------------------------------------------------%%%
%%% Lithostatic stress referenced to the top of the subducting column.
%     s_column = cumtrapz(Depth_C, g*rho);
% 
% %%% Total lithostatic stress in the column.
%     if depth(1,1) <= H
%         s_lith = s_column + interp1(Depth_C, stress_cond, depth(1,1));
%     else
% 
% %%% Vertical stress due to height of convecting ice above the top of the column.   
%         s_conv = g*rho_ice_conv.*(depth(1,1) - H);
% 
% %%% Total lithostatic stress.
%         s_lith = s_column + stress_cond(end) + s_conv;
%     end
%%%-------------------------------------------------------------------------------%%%

%%% Ice viscosity after Nimmo et al. (2003).
    eta = eta_b*exp((Q/R)*(1./T - 1/T_b));
    eta(eta > eta_max) = eta_max;

%%% Linear viscous relation for porosity change.
    dphidt = -phi.*s_lith./eta;
    dphidt(phi <= 0) = 0;

%%% Heat capacity. Johson2017 reference Kirk & Stevenson (1987)
    c_p = (1925/250)*T;                     % [J/(kg K)]

%%% Thermal conductivity. Eq. (2) in Johnson2017
    K_T = 651./T;                           % [W/(m K)]

%%% Calculate derivatives using differential operator.
    DD2 = three_point_centered_varcoeff_D2(depth(1,1), depth(N,1), N, K_T, coeff);
    dTdzz = DD2*T;

%%% Apply the top boundary condition.
    Tshell_top = min(260, T_s*(T_b/T_s)^(depth(1)/H));
    T_top = (T(2,:) + Tshell_top)/2;
    K_T_top = 651/T_top;

%%% Linear interpolation of thermal conductivity.
    dTdzz(1,:) = ((K_T_top + K_T(1))*Tshell_top - (K_T_top + 2*K_T(1)...
        + K_T(2))*T(1,:) + (K_T(1) + K_T(2))*T(2,:))/(2*dz^2);

%%% Temperature equation.
    dTdt = (1./(rho_ice.*c_p)).*dTdzz;

%%% Constant temperature at bottom boundary.
    dTdt(N,:) = 0;

%%% Assemble all of the variables.    
    dVarsdt = [dTdt; dphidt];

end %PDE

end


%%%-------------------------------------------------------------------------------%%%
%%% Create a sparsity matrix based on the Jacobian of the numerical evolution
%%% equations.
%%%-------------------------------------------------------------------------------%%%
function S = SparsityPattern(N)
    iD = (1:N);
    iSup = (1:N-1);
    I = sparse(eye(N));   

%%% Temperature evolution depends on temperature in adjacent nodes through the finite 
%%% difference formula, and does not depend on porosity.
    T_T = sparse(iD,iD,[ones(N-1,1);0]) + sparse(iSup,iSup+1,ones(N-1,1),N,N)...
        + sparse(iSup+1,iSup,[ones(N-2,1);0],N,N);
    
    Jac_T = [T_T, sparse(N,N)];

%%% Porosity evolution at each node depends on the porosity and temperature at the 
%%% same node.
    Jac_phi = [I, I];

%%% Combine all the sections to form the sparsity pattern.
    S = [Jac_T; Jac_phi];

end %SparsityPattern