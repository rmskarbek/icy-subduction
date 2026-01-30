function [R_min_factor, v_plate, SimData, Forces] = BatchRun1(R_min_factor, v_plate)

%%% This function runs a batch of simulations for different shell and slab
%%% thicknesses. Each simulation in the batch has the same plate rate, minimum plate
%%% curvature scaling factor (using the Buffett 2006 geometry) and salt structure.

%%%-------------------------------------------------------------------------------%%%
%%% INPUT
%%%-------------------------------------------------------------------------------%%%
% R_min_factor  - Scaling factor for the minimum radius of curvature in the
%                 Buffett2006 geometry. R_min = R_min_factor*H.

% v_plate       - plate convergence rate [m/yr].

%%%-------------------------------------------------------------------------------%%%
%%% Convert v_plate to m/s.
spy = 365.25*24*3600;                                       % seconds per year
v_plate_mps = v_plate/spy;                                  % [m/s]

%%% Range of shell thicknesses.
H_shell = 1e3*(10:5:70)';                                   % [m]

%%% Range of slab thicknesses for largest slab thickness.
H_70 = (500:500:H_shell(end)/2)';                           % [m]

%%% Create a table where each entry will contain a structure that contains the 
%%% simulation results and calculations.
M = numel(H_shell);
Q = numel(H_70);
sz = [Q M];

%%% Types of entries in the table.
varTypes = repmat({'struct'}, 1, M);

%%% Generate row names for slab thicknesses.
rowNames = cell(Q,1);
for i = 1:Q
    if mod(H_70(i),1e3) == 0
        rowNames{i,1} = sprintf('%d%s', 1e-3*H_70(i), 'km');
    else
        rowNames{i,1} = sprintf('%1.1f%s', 1e-3*H_70(i), 'km');
    end
end

%%% Generate column names from shell thicknesses;
varNames = cell(M,1);
for i = 1:M
    varNames{i,1} = sprintf('%d%s', 1e-3*H_shell(i), 'km');
end

%%% Table for storing the simulation data and forces calculations.
SimData = table('Size', sz, 'VariableTypes', varTypes, 'VariableNames', varNames,...
    'RowNames', rowNames);

%%% This is an annoying work around that's needed to assign structures to the table.
for i = 1:M
    SimData.(varNames{i,1}) = struct('Out', cell(Q,1), 'Bend', cell(Q,1), 'Buoyancy',...
        cell(Q,1), 'Shear', cell(Q,1));
end

%%% Generate a parameters structure will all necessary information to run a 
%%% simulation.
p = IcySubduction_Parameters;
p.Geometry.PlateRate = v_plate_mps;

%%% Allocate arrays to store the forces for each simulation.
Forces = struct('FullLength', [], 'NonSubsumedLength', []);
Forces.FullLength.F_Bend = nan(Q, M);
Forces.FullLength.F_Shear = nan(Q, M);
Forces.FullLength.F_Buoy = nan(Q, M);
Forces.NonSubsumedLength.F_Bend = nan(Q, M);
Forces.NonSubsumedLength.F_Shear = nan(Q, M);
Forces.NonSubsumedLength.F_Buoy = nan(Q, M);


%%%-------------------------------------------------------------------------------%%%
%%% Loop through the values of H_shell and H. For each pair of values, the parameters 
%%% structure needs to be updated.
%%%-------------------------------------------------------------------------------%%%
%%% Loop on H_shell.
for i = 1:M

%%% Range of slab thicknesses for current value of H_shell.
    H = (500:500:H_shell(i)/2)';                            % [m]

%%% Range of associated minimum radii of curvature.
    R_min = R_min_factor*H;                                 % [m]

%%% Loop on H.
    for j = 1:numel(H)
    
%%% Update the geometry values in the parameter structure.
        p.Geometry.ShellThick = H_shell(i);                 % [m]
        p.Geometry.SlabThick = H(j);                        % [m]
        p.Geometry.CurveRadius = R_min(j);                  % [m]

%%% Update the simulation run time, since it depends on the geometry.
        type = 'geometry';
        [~, ~, ~, ~, ~, ~, s_end] = Buffett2006(type, H(j), H_shell(i), R_min(j));
        p.Numerical.RunTime = s_end/v_plate_mps;                % [s]

%%% Update the grid spacing.
        N = p.Numerical.GridPoints;
        p.Numerical.GridSpacing = H(j)/(N-1);

%%% Update initial conditions.
        Vars_init = IcySubduction_InitialConditions(p);
        p.Numerical.InitialCond = Vars_init;

%%%-------------------------------------------------------------------------------%%%
%%% Run the simulation.
        Out = IcySubduction(p);

%%% Find the arc length location where the slab is completely subsumed.
        [ArcLengthSubsumed, k_subsumed] = SlabSubsumption(Out);
        k = numel(Out.Slab.ArcLength);

%%% Calculate the forces associated with plate bending.
        Bend = BendingForce(Out);

%%% Calculate the force due to shear resistance along the plate interface.
        mu_flag = 'temperature';
        Shear = IceShearResistance(Out, mu_flag);
        i_f = isfinite(Shear.F_Horizontal);
        S_f = Shear.F_Horizontal(i_f);
   
%%% Calculate the bouyancy force.
        Buoyancy = BuoyancyForce(Out);
    
%%%-------------------------------------------------------------------------------%%%
%%% Store the results.
        Sim = struct;
        Sim.Out = Out;
        Sim.Bend = Bend;
        Sim.Buoyancy = Buoyancy;
        Sim.Shear = Shear;
        SimData{j,i} = Sim;
        Forces.FullLength.F_Bend(j,i) = Bend.F_Fail(k);
        Forces.FullLength.F_Buoy(j,i) = Buoyancy.F_Bouy(k);
        Forces.NonSubsumedLength.F_Bend(j,i) = Bend.F_Fail(k_subsumed);
        Forces.NonSubsumedLength.F_Buoy(j,i) = Buoyancy.F_Bouy(k_subsumed);

%%% The force due to shear resistance is the same for the full length and
%%% non-subsumed length, because the shear resistance only exists along the plate
%%% interface, and the slab cannot be subsumed until it is beneath the conductive
%%% later.
        Forces.FullLength.F_Shear(j,i) = S_f(end);
        Forces.NonSubsumedLength.F_Shear(j,i) = S_f(end);
    end
end

%%%-------------------------------------------------------------------------------%%%
%%% Save the batch as a .mat file.

%%% The .mat file name is formatted as: v###_R#p#.mat. v### gives the plate rate in
%%% mm/yr. R#p# gives the R_min_factor, e.g. R_min_factor = 0.5 -> R0p5.
R_str = num2str(R_min_factor);

%%% This line to generate the .mat file name is written assuming that R_min_factor is
%%% an integer or half integer less than 10.
if isscalar(R_str)
    matName = sprintf('%s%d%s', 'v', 1e3*v_plate, ['_R' R_str(1) 'p0']);
else
    matName = sprintf('%s%d%s', 'v', 1e3*v_plate, ['_R' R_str(1) 'p' R_str(3)]);
end
save(matName, 'R_min_factor', 'v_plate', 'SimData', 'Forces', '-v7.3');