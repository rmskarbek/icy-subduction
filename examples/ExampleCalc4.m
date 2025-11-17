%%% This function will run simulations for a range of slab thicknesses and salt 
%%% contents, and a constant shell thickness. In these simulations the convective 
%%% ice and non-subducting conductive ice have no salt; alternatively, the salt 
%%% content of the slab can be interpreted as the difference between the slab and the 
%%% other ice.

%%% The function creates two plots for the set of simulations: one where the forces
%%% are computed along the entire length of the slab (i.e. up until the slab hits the
%%% ice-ocean interface), and a second plot where the forces are computed only along
%%% the length of the slab that is not subsumed.

%%%-------------------------------------------------------------------------------%%%

%%% Range of slab thicknesses.
% H = 1e3*[0.5; (1:8)'];                                         % [m]
H = 1e3*(0.5:0.5:4)';
% H = 1e3*(0.5:0.5:6)';

%%% Range of associated minimum radii of curvature.
R_min = 2*H;                                            % [m]

%%% Range of slab salt content.
% Salt_slab = [0; 0.05; 0.1; 0.15; 0.2];
Salt_slab = (0:0.025:0.2)';

%%% Create a table where each entry will contain a structure that contains the 
%%% simulation results and calculations.
M = numel(H);
Q = numel(Salt_slab);
sz = [Q M];

%%% Types of entries in the table.
varTypes = repmat({'struct'}, 1, M);

%%% Generate row names for salt contents.
rowNames = cell(Q,1);
for i = 1:Q
    rowNames{i,1} = sprintf('%d%s', 1e2*Salt_slab(i), '%');
end

%%% Generate column names from slab thicknesses;
varNames = cell(M,1);
for i = 1:M
    varNames{i,1} = sprintf('%d%s', 1e-3*H(i), 'km');
end

T = table('Size', sz, 'VariableTypes', varTypes, 'VariableNames', varNames,...
    'RowNames', rowNames);

%%% This is an annoying work around that's needed to assign structures to the table.
for i = 1:M
    T.(varNames{i,1}) = struct('Out', cell(Q,1), 'Bend', cell(Q,1), 'Buoyancy',...
        cell(Q,1), 'Shear', cell(Q,1));
end

%%% Generate a parameters structure will all necessary information to run a 
%%% simulation.
p = IcySubduction_Parameters;
% p.Geometry.PlateRate = v_plate;

%%% Allocate arrays to store the forces for each simulation.
Forces = struct('FullLength', [], 'NonSubsumedLength', []);
Forces.FullLength.F_Bend = nan(Q, M);
Forces.FullLength.F_Shear = nan(Q, M);
Forces.FullLength.F_Buoy = nan(Q, M);
Forces.NonSubsumedLength.F_Bend = nan(Q, M);
Forces.NonSubsumedLength.F_Shear = nan(Q, M);
Forces.NonSubsumedLength.F_Buoy = nan(Q, M);

%%%-------------------------------------------------------------------------------%%%
%%% Loop through the values of H and F_slab. For each value, the parameters structure 
%%% needs to be updated.
%%%-------------------------------------------------------------------------------%%%
for i = 1:M
    for j = 1:Q
    
%%% Update the slab thickness and minimum radius of curvature.
        p.Geometry.SlabThick = H(i);                        % [m]
        p.Geometry.CurveRadius = R_min(i);                  % [m]

%%% Update the simulation run time, since it depends on the geometry.
        type = 'geometry';
        H_shell = p.Geometry.ShellThick;                    % [m]
        v_plate = p.Geometry.PlateRate;                     % [m/s]

        [~, ~, ~, ~, ~, ~, s_end] = Buffett2006(type, H(i), H_shell, R_min(i));
        p.Numerical.RunTime = s_end/v_plate;                % [s]

%%% Update the grid spacing.
        N = p.Numerical.GridPoints;
        p.Numerical.GridSpacing = H(i)/(N-1);

%%% Update the salt content in the slab.
        p.Porosity.SaltSlab = Salt_slab(j);

%%% Update initiatl conditions.
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
        % mu_flag = 'constant';
        mu_flag = 'temperature';
        Shear = IceShearResistance(Out, mu_flag);
        i_f = isfinite(Shear.F_Horizontal);
        S_f = Shear.F_Horizontal(i_f);
    
%%% Use these lines if you want to use only the frictional resistance on the plate
%%% interface.
        % i_f = isfinite(Shear.F_Frict);
        % S_f = Shear.F_Frict(i_f);

%%% Calculate the bouyancy force.
        Buoyancy = BuoyancyForce(Out);
    
%%%-------------------------------------------------------------------------------%%%
%%% Store the results.
        Sim = struct;
        Sim.Out = Out;
        Sim.Bend = Bend;
        Sim.Buoyancy = Buoyancy;
        Sim.Shear = Shear;
        T{j,i} = Sim;
        Forces.FullLength.F_Bend(j,i) = Bend.F_Fail(k);
        Forces.FullLength.F_Buoy(j,i) = Buoyancy.F_Bouy(k);
        Forces.NonSubsumedLength.F_Bend(j,i) = Bend.F_Fail(k_subsumed);
        Forces.NonSubsumedLength.F_Buoy(j,i) = Buoyancy.F_Bouy(k_subsumed);

%%% This is meant to ensure that the shear resistance is only computed along the
%%% non-subsumed length of the slab. But it is unecessary, since the slab cannot be
%%% subsumed along the plate interface.
        % F_Shear(j,i) = S_f(min(k, numel(S_f)));
        Forces.FullLength.F_Shear(j,i) = S_f(end);
        Forces.NonSubsumedLength.F_Shear(j,i) = S_f(end);
    end
end

%%%-------------------------------------------------------------------------------%%%
%%% Plot the results for both the full length of the slab and the non-subsumed
%%% length.
for i = 1:numel(fieldnames(Forces))
    if i == 1
        F_tot = Forces.FullLength.F_Buoy./(Forces.FullLength.F_Bend...
            + Forces.FullLength.F_Shear);
    else
        F_tot = Forces.NonSubsumedLength.F_Buoy./(Forces.NonSubsumedLength.F_Bend...
            + Forces.NonSubsumedLength.F_Shear);
    end

    [X, Y] = meshgrid(1e-3*H, Salt_slab);
    Slabs = reshape(X, [numel(F_tot) 1]);
    Salts = reshape(Y, [numel(F_tot) 1]);
    ForcesVec = reshape(F_tot, [numel(F_tot) 1]);
    
    [H_grid, Salt_grid] = meshgrid(1e-3*linspace(H(1), H(end), 30),...
        linspace(Salt_slab(1), Salt_slab(end), 30));
    Forces_grid = griddata(Slabs, Salts, ForcesVec, H_grid, Salt_grid);
    
    font = 'Palatino Linotype';
    f = figure('DefaultTextFontName',font,'DefaultAxesFontName',font,...
        'DefaultAxesFontSize', 16);
    % theme(f, "dark")
    theme(f, "light")
    
    % [x1, ~] = contour(H_grid, Salt_grid, Forces_grid, [1 1], 'w', 'LineWidth', 2);
    [x1, ~] = contour(H_grid, Salt_grid, Forces_grid, [1 1], 'k', 'LineWidth', 2);
    hold on
    scatter(Slabs, Salts, 200, ForcesVec, 'filled')
    hold off
    
    ax = gca;
    ax.XLim = [0 10.5];
    ax.YLim = [-.01 0.21];
    ax.TickLabelInterpreter = 'latex';
    xlabel('Slab Thickness (km)', 'Interpreter', 'latex')
    ylabel('Slab Salt Content', 'Interpreter', 'latex')
    
    c = colorbar;
    c.TickLabelInterpreter = 'latex';
    c.Title.Interpreter = 'latex';
    c.Title.String = '$F_{drive}/F_{resist}$';
    c.Title.FontSize = 18;
    
    % cmap = colormap(ax, 'winter');
    cmap = crameri('roma');
    if i == 1
        y = (1 - min(min(F_tot)))/(max(max(F_tot)) - min(min(F_tot)));
        % y = 1/(max(max(F_tot)) - min(min(F_tot)));
        x = round((256*y - 128)/(y-1));
        cmap2 = cmap(x:end,:);
        cmap2(1:round(y*length(cmap2)),:) = cmap(round(linspace(1, 256/2,...
            numel(1:round(y*length(cmap2))))'),:);
        ForcesVec1 = ForcesVec;
        colormap(ax, cmap2);
    else
        i_map3 = round((max(ForcesVec)/max(ForcesVec1))*size(cmap2,1));
        cmap3 = cmap2(1:i_map3,:);
        colormap(ax, cmap3);
    end        
end