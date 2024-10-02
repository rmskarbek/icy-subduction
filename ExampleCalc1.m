%%% This script uses temperature profiles and shell thicknesses from Walker and Rhoden
%%% (2022) to calculate the resistance to sliding on a inclined fault that goes through
%%% the ice shell.

%%% First, get the temperature data as a function of depth for an eccentricity of 0.01, as
%%% well as the shell thicknesses. This will also produce a plot of the temperature
%%% profiles.
[Depth, Shells, T] = ReadTemperatureData;

%%% Second, define the range of dip angles for the calculation. Must be in degrees.
DipAngle = (10:1:45)';

%%% Force per unit distance to cause sliding on the fault plane. There will be one vector
%%% for each value of shell thickness.
Resistance = nan(numel(DipAngle), numel(Shells));

%%% Flags for the type of calculate and how to handle the friction coefficient.
calc = 'fault';
mu_flag = 'temperature';

%%% Run the calculation.
for j = 1:numel(Shells)
    for i = 1:numel(DipAngle)

    %%% Get the temperature profile and remove nan entries.
        TempK = T{:,j};
        i_finite = isfinite(TempK);
        TempK = TempK(i_finite);

    %%% Get the corresponding depth.
        Depth_shell = Depth(i_finite);

    %%% 9.19.2024 - Should also store the shear resistance values in a table or something.
        [~, ~, S, F] = IceShearResistance(Depth_shell, TempK, calc, mu_flag, DipAngle(i));       
        Resistance(i, j) = F;
           
    end
end

%%% Plot all the resistances.
figure;
hold on;
box on
for i = 1:numel(Shells)
    plot(DipAngle, Resistance(:,i), 'LineWidth', 2)
end
hold off

ax = gca;
ax.YScale = 'log';

xlabel('Fault Dip Angle (deg)')
ylabel('Shear Resistance (N/m)')
title('Thinner shell means less resistance')

lgd = legend(cellstr(num2str(Shells)), 'Location', 'eastoutside');
title(lgd, 'Shell Thickness (km)')