function [Depth, Shells, T] = ReadTemperatureData

%%% This script will read temperature profile data from Rhoden and Walker (2022).
filename = 'PSIE_TempProfiles_Sept2024.xlsx';

%%% Import the data from a specified sheet.
sheet = 'ecc = 0.0100';
T1 = readtable(filename, 'Sheet', sheet);

%%% Create a new table that is just the temperature values.
T = T1(4:end, 2:end);

%%% Label each row in the temperature table with the corresponding depth.
Depth = T1{4:end, 1};                                       % [km]
T.Properties.RowNames = cellstr(num2str(Depth));

%%% Label each column in the temperature table with the corresponding shell thickness.
Shells = T1{1, 2:end}';
ShellNames = cellstr(num2str(Shells))';
T.Properties.VariableNames = ShellNames;

%%% Delete T1 and remove unused rows from T and Depth.
clear T1
k = isfinite(T{:,end});
N = numel(T{:,end}(k));
T = T(1:N,:);
Depth = Depth(1:N);

%%% Plot all the temperature profiles.
figure;
hold on;
box on
for i = 1:numel(T.Properties.VariableNames)
    plot(T{:,i}, 1e-3*Depth, 'LineWidth', 2)
end
ax = gca;
ax.YDir = 'reverse';

xlabel('Temperature (K)')
ylabel('Depth (km)')
title('Temperature Profiles from Rhoden and Walker (2022); ecc = 0.01')
lgd = legend(ShellNames, 'Location', 'southwest');
title(lgd, 'Shell Thickness (km)')