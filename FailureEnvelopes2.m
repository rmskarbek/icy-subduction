function [Depth, DipAngle, Resistance, S_30] = FailureEnvelopes2

%%% Dip angle of subduction fault plane [radians]. See Fig. 1 in Mueller & Phillips.
DipAngle = (pi/180)*(10:0.1:45)';
[~, i_30] = min(abs((180/pi)*DipAngle - 30));

%%% Set up a figure for plotting the results.
font = 'DejaVu Serif';
figure('DefaultTextFontName', font, 'DefaultAxesFontName', font,...
    'DefaultAxesFontSize',  14)
set(gcf, 'Color', 'w')

%%% Force per unit distance to cause sliding on the fault plane.
Resistance = nan(size(DipAngle));

%%% Flags for how to handle the friction coefficient and line styles for each case.
calc = 'fault';
Mu_Flags = {'temperature'; 'constant'};
LineStyle = {'-'; '-.'};
for j = 1:numel(Mu_Flags)
    for i = 1:numel(DipAngle)
        [Xi, Depth, TempK, S, F] = IceShearResistance(DipAngle(i), calc, Mu_Flags{j});
        Resistance(i) = F;
    
        if i == i_30
        %%% Plot the shear resistance against the depth on the fault plane.
            subplot(1, 2, 1);
            hold on
            S_30 = S;
            plot(S, Depth, LineStyle{j}, 'LineWidth', 2)
            hold off
        end
    end
end

%%% Labels and such for the shear resistance plot.
box on
grid on
ax = gca;
ax.YDir = 'reverse';
ax.FontSize = 14;
title('Europa Shear Resistance, 30 deg dip')
xlabel('Shear Resistance (MPa)', 'FontSize', 16)
ylabel('Depth (km)', 'FontSize', 16)
legend('Temperature dependent', '\mu = 0.55, Howell & Pappalardo (2019)');

%%% Plot the temperature profile next to the shear resistance profile.
subplot(1, 2, 2)
plot(TempK, Depth, 'k-', 'LineWidth', 2)
grid on
ax = gca;
ax.YDir = 'reverse';
ax.FontSize = 14;
% title('Europa Shear Resistance, 30 deg dip')
xlabel('Temperature (k)', 'FontSize', 16)
ylabel('Depth (km)', 'FontSize', 16)

%%% Plot the resistance against the dip angle.
% figure
% plot((180/pi)*DipAngle, 1e-9*Resistance, 'LineWidth', 2)
% grid on
% ax = gca;
% ax.FontSize = 14;
% ax.XLim = [0 50];
% % ax.YLim = [1e13 16e13];
% title('Europa, Shear Resitance Force')
% xlabel('Fault Plane Dip (deg)', 'FontSize', 16)
% ylabel('Resistance (GN/m)', 'FontSize', 16)

%%% Output dip angle in degrees.
DipAngle = (180/pi)*DipAngle;

end