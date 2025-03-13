function Out = IcySubduction(p)
%%%------------------------------------------------------------------------------------%%%
%%%------------------------------------------------------------------------------------%%%
%%% Notes.
%%% 1. Account for density variations of ice above the slab in sigma_L? In
%%% IceFailureEnvelope.m and IceShearResistance.m. 
 
%%% 2. Also reconcile the temperature dependent densities in IceFailureEnvelope.m and 
%%% IceShearResistance.m with the porosity output.

%%% 3. Move numeric parameters (grid points, run time) to the parameters structure in
%%% IcySubduction_Parameters.m.

%%% 4. Stop the simulation when min(Temp) = BasalTemp & max(Phi) = 0.

%%% 5. Make separate initial condtions code that takes parameters structure.

%%% 6. Update the Johnson geometry calculations.

%%%------------------------------------------------------------------------------------%%%
%%% This code solves the dimensional equations for temperature and porosity evolution in a
%%% subducting 1-D column of ice. The code is essentially the same model used by
%%% Johnson et al. (2017) and Howell & Papallardo (2019), although there are some
%%% important differences:

%%% 1. The code includes an option to use the slab geometry defined by Buffett (2006).
%%%    Whereas J2017 and H2019 both used a circular arc slab geometry.

%%% 2. This code employs the correct finite difference formula for dealing with
%%%    non-constant coefficients in the diffusive term of the temperature equation. Both
%%%    J2017 and H2019 use an incorrect formula.

%%%------------------------------------------------------------------------------------%%%
%%% INPUT
% p   - A matlab structure containting parameter values for the simulation. Generate with:
%       p = IcySubduction_Parameters.m

%%% OUTPUT
% Out - A structure that contains the simulations output, the slab geometry, and also the
%       same parameters structure that is input.


%%%------------------------------------------------------------------------------------%%%
%%% Constants.
%%%------------------------------------------------------------------------------------%%%
spy = 365.25*24*3600;                   % seconds per year
g = p.Constants.Gravity;                % Europa gravity [m/s^2]
R = p.Constants.GasConst;               % gas constant [J/(mol K)]
rho_ice = p.Constants.IceDensity;       % reference ice density [kg/m^3]

%%%------------------------------------------------------------------------------------%%%
%%% Geometry
%%%------------------------------------------------------------------------------------%%%
H = p.Geometry.SlabThick;               % thickness of conductive slab [m]
H_shell = p.Geometry.ShellThick;        % total thickness of ice shell [m]
GeoFlag = p.Geometry.GeoFlag;           % a flag that determines the slab geometry.
R_min = p.Geometry.CurveRadius;         % minimum radius of plate curvature [m]
v_plate = p.Geometry.PlateRate;         % plate convergence rate [m/s]
type = 'location';                      % flag for Buffett2006.m

%%%------------------------------------------------------------------------------------%%%
%%% Temperature
%%%------------------------------------------------------------------------------------%%%
T_b = p.Temperature.BasalTemp;          % basal temp [K] 
T_s = p.Temperature.SurfaceTemp;        % surface temp [K], after Nimmo et al. (2003)

%%%------------------------------------------------------------------------------------%%%
%%% Viscosity
%%%------------------------------------------------------------------------------------%%%
Q = p.Viscosity.ActivEnergy;            % viscosity activation energy [J/mol]
eta_b = p.Viscosity.BasalVisc;          % reference viscosity at base of ice shell [Pa*s]
eta_max = p.Viscosity.MaxVisc;          % maximum viscosity from Howell2019 [Pa*s]

%%%------------------------------------------------------------------------------------%%%
%%% Porosity
%%%------------------------------------------------------------------------------------%%%
t_por = p.Porosity.CompTime;            % initial porosity compaction time [s]
phi_0 = p.Porosity.PorosityRef;         % porosity reference


%%%------------------------------------------------------------------------------------%%%
%%% Numerical
%%%------------------------------------------------------------------------------------%%%
N = p.Numerical.GridPoints;             % Number of grid points
dz = p.Numerical.GridSpacing;           % Numerical grid spacing [m]

% RunTime = p.Numerical.RunTime;          % run time of the simulation [s]
RunTime = 550e3*spy;                    % run time of the simulation [s]
time = linspace(0, RunTime, 1e3);       % vector of time values for solution output.


%%%------------------------------------------------------------------------------------%%%
%%% Initial conditions.
%%%------------------------------------------------------------------------------------%%%
Depth_init = linspace(0, H, N)';        % depth referenced to the top of the column [m]

%%% Initial lithostatic stress.
sigma_lith = g*rho_ice*Depth_init;                           % lithostatice stress [Pa]

%%% Initial temperature profile in the subducting slab. Eq. (3) in Johnson2017
T_init = T_s*(T_b/T_s).^(Depth_init/H);                     % initial temperature [K]

%%% Initial viscosity profile in the subducting slab. Eq. (5) in Johnson2017
eta_init = eta_b*exp((Q/R)*(1./T_init - 1/T_b));            % initial viscosity [Pa s]

%%% Initial porosity profile in the subducting slab. Eq. (4) in Johnson2017
phi_init = phi_0*exp(-t_por*sigma_lith./eta_init);

%%% Construct the intial conditions vector.
Vars_init = [T_init; phi_init];


%%%------------------------------------------------------------------------------------%%%
%%% Differential operators and sparsity pattern.
%%%------------------------------------------------------------------------------------%%%
%%% Differentiation matrices for finite differences.
%%% Set v = -1 to indicate that waves are moving towards shallower depths 
%%% (i.e. in negative z direction).
    % v = -1;
%%% Three point upwind difference for first derivative.
    % D1 = three_point_upwind_uni_D1(Depth(1,1), Depth(N,1), N, v);

%%% Three point centered difference for second derivative.
    % DD2 = three_point_centered_uni_D2(Depth_init(1,1), Depth_init(N,1), N);

%%% Three point centered difference for second derivative, variabl coefficients.
    % DD2 = three_point_centered_varcoeff_D2(Depth_init(1,1), Depth_init(N,1), N, K_T_init);

%%% Generate the sparsity pattern for amazing speed up!
    % S = SparsityPattern(N);
    % options = odeset('Vectorized','on','JPattern',S,'RelTol',1e-10,'AbsTol',...
    %     [1e-12*ones(N,1); 1e-12*ones(N,1); 1e-12*ones(N,1); 1e-12*ones(N,1)]);
%%% Solve the system.    
    % [Time, Vars_f] = ode15s(@PDE, Time, Vars_i, options);

    
%%%------------------------------------------------------------------------------------%%%
%%% Conduct the simulation and deal with the output.
%%%------------------------------------------------------------------------------------%%%
[Time, Vars_f] = ode15s(@PDE, time, Vars_init);
% [Time, Vars_f] = ode89(@PDE, time, Vars_init);

Temp = Vars_f(:, 1:N)';
Phi = Vars_f(:, N+1:2*N)';
ArcLength = Time*v_plate;                                                       % [m]

%%% Evalute the geometry for the output time.
switch GeoFlag
    case 'Buffett'
        [Z_Top, Z_Center, Z_Bottom, Theta, K, dKds] = Buffett2006(type, H, H_shell, R_min,...
            ArcLength);

    case 'Johnson'
        [Z_Top(i), Z_Center(i), Z_Bottom(i)] = JohnsonGeometry(Time(i));
end

%%% Get the depth of the slab column at each output time.
Depth = nan(N, numel(Time));
for i = 1:numel(Time)
    Depth(:, i) = linspace(imag(Z_Top(i,1)), imag(Z_Bottom(i,1)), N)';
end

%%% Convert time to years.
Time = Time/spy;

%%% Evaluate the viscosity in the slab.
Eta = eta_b*exp((Q/R)*(1./Temp - 1/T_b));
Eta(Eta > eta_max) = eta_max;

%%% Create output structure.
Out = struct('TimeYrs', Time, 'Temperature', Temp, 'Porosity', Phi, 'Viscosity', Eta,...
    'Slab', [], 'p', p);
Out.Slab.ArcLength = ArcLength;
Out.Slab.Depth = Depth;
Out.Slab.CenterLine = Z_Center;
Out.Slab.Top = Z_Top;
Out.Slab.Bottom = Z_Bottom;
Out.Slab.Dip = Theta;
Out.Slab.Curvature = K;
Out.Slab.dKds = dKds;

%%%------------------------------------------------------------------------------------%%%
%%%------------------------------------------------------------------------------------%%%
function dVarsdt = PDE(time, Vars, ~)
    T = Vars(1:N,:);
    phi = Vars(N+1:2*N,:);

%%% Get the current burial depth of the numerical grid points.
    switch GeoFlag
        case 'Buffett'
            s = v_plate*time;
            [z_top, ~, z_bottom] = Buffett2006(type, H, H_shell, R_min, s);

        case 'Johnson'
            [z_top, ~, z_bottom] = JohnsonGeometry(time);
    end
    depth = linspace(imag(z_top), imag(z_bottom), N)';

%%% Lithostatic stress.
    s_lith = g*rho_ice*depth;

%%% Upper boundary.
    % Tshell_top = min(260, T_s*(T_b/T_s)^(Depth(1)/h_cond));
    % T(1) = Tshell_top;
    % T_top = (T(2) + Tshell_top)/2;
    % T(1) = T_top;

%%% Ice viscosity after Nimmo et al. (2003).
    eta = eta_b*exp((Q/R)*(1./T - 1/T_b));
    eta(eta > eta_max) = eta_max;

%%% Linear viscous relation for porosity change dphi = -phi*p_eff/eta.
    % p_eff = repmat(p_lith + dsigmadt*time,1,size(phi,2)) - p_tot;
    dphidt = -phi.*s_lith./eta;
    dphidt(phi <= 0) = 0;

%%% Temperature equation. Update the thermal conducticity and differential operator.
%%% Thermal conductivity. Eq. (2) in Johnson2017
    K_T = 651./T;                           % [W/(m K)]
    DD2 = three_point_centered_varcoeff_D2(depth(1,1), depth(N,1), N, K_T);

%%% Calculate derivatives using differential operator.
    % dTdz = D1*T;
    dTdzz = DD2*T;

%%% Apply the top boundary condition.
    Tshell_top = min(260, T_s*(T_b/T_s)^(depth(1)/H));
    % T_top = Tshell_top;
    T_top = (T(2,:) + Tshell_top)/2;
    K_T_top = 651/T_top;

    dTdzz(1,:) = ((K_T_top + K_T(1))*Tshell_top - (K_T_top + 2*K_T(1) + K_T(2))*T(1,:)...
        + (K_T(1) + K_T(2))*T(2,:))/(2*dz^2);

    % dTdzz(1,:) = (Tshell_top - 2*T(1,:) + T(2,:))/(2*dz^2);

%%% Heat capacity. Johson2017 reference Kirk & Stevenson (1987)
    c_p = (1925/250)*T;                     % [J/(kg K)]

%%% Only heat conduction.
    % dTdt = (K_T./(rho_ice*c_p)).*dTdzz;
    dTdt = (1./(rho_ice*c_p)).*dTdzz;

%%% Set zero change at top boundary, but really this should be changed to use a ghost
%%% node.
    % dTdt(1,:) = 0;

%%% Constant temperature at bottom boundary.
    dTdt(N,:) = 0;

%%% Assemble all of the variables.    
    dVarsdt = [dTdt; dphidt];

end %PDE


%%%------------------------------------------------------------------------------------%%%
%%%------------------------------------------------------------------------------------%%%
function [z_top, z_center, z_bottom] = JohnsonGeometry(time)

    DipAngle = 45*pi/180;               % Final slab dip angle [radians]
%%% Current arc length, tangent angle, and coordinates along the center line.
    s = v_plate*time;

    if time <= R_min*DipAngle/v_plate
        theta = s/R_min;
        z_center = 1i*R_min*(1 - exp(1i*s/R_min)) + 1i*H/2;
    else
        theta = DipAngle;
        z_center = 1i*R_min*(1 - exp(1i*DipAngle)) + exp(1i*DipAngle)*(s - R_min*DipAngle)...
            + 1i*H/2;
    end

%%% Finally, compute the depth of the slab upper surface by projecting a distance h_cond/2
%%% along the angle normal to the center line.
%%% Normal angle (psi) along the center line.
    psi = pi/2 - theta;

%%% Slope of normal line.
    m = -tan(psi);

%%% y-intercept of normal line.
    % b = imag(z_center) - m.*real(z_center);
    x_center = real(z_center);
    y_center = imag(z_center);
    b = y_center - m.*x_center;    

%%% Coordinates of the top and bottom surfaces.
    A = (1 + m.^2);
    % B = 2*(m.*(b - imag(z_center)) - real(z_center));
    % C = real(z_center).^2 + (b - imag(z_center)).^2 - (H/2)^2;
    B = 2*(m.*(b - y_center) - x_center);
    C = x_center.^2 + (b - y_center).^2 - (H/2)^2;
    
    x_top = (1./(2*A)).*(-B + sqrt(B.^2 - 4*A.*C));
    x_bottom = (1./(2*A)).*(-B - sqrt(B.^2 - 4*A.*C));
    
    y_top = m.*x_top + b;
    y_bottom = m.*x_bottom + b;
    
    z_top = x_top + 1i*y_top;
    z_bottom = x_bottom + 1i*y_bottom;

end

end

%%%------------------------------------------------------------------------------------%%%
%%%------------------------------------------------------------------------------------%%%
function [Pressure_i,Pressure_0,Hydrostatic,Temp,dTempdt,dsigmadt,dp_lithdz]...
    = Initial(phi_i,N,Depth,Depth_0,rho_s,rho_f)
%%% This function computes the initial conditions in the column based on
%%% the starting depth of the simulation and the total depth. And also
%%% computes some parameter values related to ongoing subduction.

%%% INPUT
% phi_i                         Initial porosity in the column [-]
% N                             Number of nodes in the column
% Depth                         Depth referenced to the top of the column [m]
% Depth_0                       Depth of the top of the column referenced
%                               to the surface [m]
% rho_f                         fluid density [kg/m^3]
% rho_s                         rock density [kg/m^3]

%%% OUTPUT
% Pressure_i                    Excess fluid pressure in column [Pa]
% Pressure_0                    Total fluid pressure in column [Pa]
% Hydrostatic                   Hydrostatic fluid pressure in column [Pa]
% Temp                          Temperature in column [K]
% dTempdt                       Rate of temperature increase at the base of
%                               the column due to subduction [K/s]
% dsigmadt                      Rate of increase in the overburden due to
%                               subduction [Pa/s]
% d_plithdz                     Lithostatic overpressure gradient [Pa/m]


%%% Some examples:
%%% Initial pressure and temperature profiles. 
%%% For Depth = 1000 m:
%%% 0 nodes dehydrating: 32.19355e3
%%% 1 node dehydrating: 32.1962e3
%%% About 12 m dehydrating: 32.24e3
%%% About 50 nodes dehydrating: 32.4e3

%%% For Depth  = 2000 m
%%% About 50 m dehydrating: 28.5e3
%%% About 12 m dehydrating: 28.345e3, gives same DeltaG_i profile as
%%% 32.24e3 for 1000 m column.

%%% For Depth = 3000
%%% About 12 m dehydrating: 24.383e3, gives same DeltaG_i profile as
%%% 32.24e3 for 1000 m column.

%%% Plate interface for distances far from the trench. The coefficients of 
%%% the plate interface polynomial are found from a linear fit to the plate
%%% interface depth for depths ~20-30 km
    PlateInt = [0.2709 -9.8967e3]; %(meters)%
    Distance_0 = (Depth_0 - PlateInt(2))/PlateInt(1);
    
%%% Geometrical Parameters for Temperature and LoadingConstant calculation.
    betaT = atan(PlateInt(1));
    v_plate = 0.037/(365*24*3600); % plate convergence velocity [m/s]
    v_x = v_plate*cos(betaT);
    StressFit = [8693.11664571161 -365056544.708814];
    dsigmadt = StressFit(1)*cos(betaT)*v_plate;
   
%%% Temperature field estimated from Peacock 2009. TempDepth is the linear
%%% increase with subduction depth.
    TempDepth = [200/30000 (400*50000-600*20000)/30000];
    dTdz = 500/25000;
    dTempdt = TempDepth(1)*PlateInt(1)*v_x;
    Temp = TempDepth(1).*(PlateInt(1).*Distance_0+PlateInt(2))+TempDepth(2)...
        + dTdz.*Depth + 273;
%%% Calculate the stress fit and initial overburden. The stress coefficients
%%% are found from a linear fit the stress estimated at the plate interface
%%% from the Trehu velocity model, for the same range as the plate interface
%%% polynomial.
    StressFit = [8693.11664571161 -365056544.708814]; %(Pa%
    Overburden = polyval(StressFit,Distance_0);      
%%% Set the initial porosity.        
    Porosity = phi_i.*ones(N,1);
    
%%% Calculate the lithostatic stress in the column due to the self-weight,
%%% assuming an initial exess pore pressure, Pressure_0, in the column.
    Lithostatic = zeros(N,1);
    BulkDensity = rho_s - (rho_s-rho_f).*Porosity; 
    Lithostatic(1,1) = 9.81*(Depth(1))*BulkDensity(1);
        for i = 2:N
            Lithostatic(i,1) = 9.81*(Depth(i) - Depth(i-1))*BulkDensity(i)...
                + Lithostatic(i-1);                
        end    
%%% Hydrostatic pressure with reference datum set to the top of the column.    
    Hydrostatic = 9.81*rho_f*(Depth+Depth_0);
%%% This is the total pressure, the term proportional to the overburden is
%%% the initial excess pressure. For lithostatic pore pressure, the EXCESS
%%% pore pressure P_e = Overburden + Lithostatic - Hydrostatic; thus the 
%%% TOTAL pore pressure is P_t = Overburden + Lithostatic, and the
%%% effective stress will then be zero everywhere.
%%% Pressure_0 is the TOTAL initial pore pressure.
    Pressure_0 = (Overburden).*ones(N,1) + Lithostatic;
    Pressure_i = Pressure_0 - Hydrostatic;
%%% Calculate the lithostatic overpressure gradient.    
    dp_lithdz = 9.81*(BulkDensity - rho_f);
    dp_lithdz = mean(dp_lithdz);
end %Initial

%%%------------------------------------------------------------------------------------%%%
%%%------------------------------------------------------------------------------------%%%
function S = SparsityPattern(N)
%%% Numerical set up. Create the JPattern sparsity matrix based on the
%%% finite-difference equations.
    iD = (1:N);
    iSup = (1:N-1);
    iSupSup = (1:N-2);
%%% Pressure section.
    P_p = sparse(iD,iD,ones(N,1)) + sparse(iSup,iSup+1,ones(N-1,1),N,N)...
        + sparse(iSup+1,iSup,ones(N-1,1),N,N)...
        + sparse(iSupSup,iSupSup+2,ones(N-2,1),N,N);
    P_p(N,N-2) = 1;
    P_p(N,N-3) = 1;
    P_T = sparse(iD,iD,ones(N,1));
    P_phi = sparse(iD,iD,ones(N,1),N,N) + sparse(iSup,iSup+1,ones(N-1,1),N,N)...
        + sparse(iSupSup,iSupSup+2,ones(N-2,1),N,N);
    P_phi(N-1,N-2) = 1;
    P_phi(N,N-1) = 1;
    P_phi(N,N-2) = 1;
    P_m = P_T;
    Jac_P = [P_p P_T P_phi P_m];
%%% Temperature equation. 
    T_p = sparse(iD,iD,[ones(N-1,1);0]) + sparse(iSup,iSup+1,ones(N-1,1),N,N)...
        + sparse(iSup+1,iSup,[ones(N-2,1);0],N,N);
    T_T = T_p;
    T_phi = P_phi;
    T_m = sparse(iD,iD,[ones(N-1,1);0]);
    Jac_T = [T_p T_T T_phi T_m];
%%% Porosity section is equal to the pressure section.
%%% Dehydration section.
    Jac_m = [P_T P_T sparse(N,N) P_T];
%%% Combine all the sections to form the sparsity pattern.
    S = [Jac_P; Jac_T; Jac_P; Jac_m];
end %SparsityPattern
