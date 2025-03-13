function GeometryPlot(p)

%%%------------------------------------------------------------------------------------%%%
%%% This function plots the model geometry with the slab, ice shell and coordinate axes.
%%% It also annotates the figure to the show the temperature boundary conditions.

%%% NOTE: This code is optimized for H = 5 km and H_shell = 25 km, but should still work
%%% for other thicknesses.
%%%------------------------------------------------------------------------------------%%%

%%% Get the geometry of the slab.
type = 'geometry';
H = p.Geometry.SlabThick;
H_shell = p.Geometry.ShellThick;
R_min = p.Geometry.CurveRadius;
GeoFlag = 'Buffett';                % a flag that determines the slab geometry.

switch GeoFlag
    case 'Buffett'
        [w_top, w_center, w_bottom, theta, ~, ~, ~] = Buffett2006(type, H, H_shell, R_min);

    case 'CircularArc'
        3;
end

%%% Convert units to kilometers.
Slab = 1e-3*[w_bottom, w_center, w_top];
H = 1e-3*H;
H_shell = 1e-3*H_shell;

%%% Set up the figure.
font = 'Palatino Linotype';
fig = figure('DefaultTextFontName',font,'DefaultAxesFontName',font,...
    'DefaultAxesFontSize', 16);
set(gcf,'Color','w')

xf = ceil(max(H_shell, max(real(Slab(:,3)))));
x1 = -10;
% x1 = -xf;
x2 = xf;

%%% Plot the convective layer.
area([x1 x2], H_shell*[1 1], FaceColor = "#0072BD");
hold on

%%% Plot the conductive later.
area([x1 x2], H*[1 1], FaceColor = 'b');

%%% Plot the slab.
fill(real([x1; Slab(:,3); flip(Slab(:,1)); x1]), imag([0; Slab(:,3); flip(Slab(:,1));...
    1i*H]), 'b')

%%%------------------------------------------------------------------------------------%%%
%%% Plot cartesian coordinate axes.
Arrow_X = annotation('arrow') ;
Arrow_X.Parent = gca;
Arrow_X.Position = [0, 0, 10, 0];
Arrow_X.LineWidth = 2;

Arrow_Y = annotation('arrow') ;
Arrow_Y.Parent = gca;
Arrow_Y.Position = [0, 0, 0, 10];
Arrow_Y.LineWidth = 2;
Arrow_Y.HeadLength = -10;

text(10, -0.5, '$x$', 'FontSize', 20, 'Interpreter', 'latex')
text(0.5, 10, '$y$', 'FontSize', 20, 'Interpreter', 'latex')

%%%------------------------------------------------------------------------------------%%%
%%% Plot the center line and the trench location.
plot(real(Slab(:,2)), imag(Slab(:,2)), 'k--', LineWidth = 2)
plot(real([Slab(1,1) Slab(1,3)]), imag([Slab(1,1) Slab(1,3)]), 'k', LineWidth = 2)
plot(real([Slab(end,1) Slab(end,3)]), imag([Slab(end,1) Slab(end,3)]), 'k', LineWidth = 2)

%%% Plot the plate interface.
[~, i_int] = min(abs(imag(Slab(:,3)) - H));
i_int = i_int + 1;
plot(real(Slab(1:i_int,3)), imag(Slab(1:i_int,3)), 'r--', LineWidth = 2)

% %%% Plot some locations of the column.
% plot(real(Slab(1, [1 3])), imag(Slab(1, [1 3])), 'k', LineWidth=3)
ic2 = round(0.3*numel(theta));
% plot(real(Slab(ic2, [1 3])), imag(Slab(ic2, [1 3])), 'k', LineWidth=3)
% ic3 = round(0.9*numel(theta));
% plot(real(Slab(ic3, [1 3])), imag(Slab(ic3, [1 3])), 'k', LineWidth=3)
%

%%%------------------------------------------------------------------------------------%%%
%%% Tangent angle (psi) along the center line.
psi = pi/2 - theta;
%%% Slope of normal lines.
m = -tan(psi);
%%% y-intercept of normal lines.
y_int = imag(Slab(:,2)) - m.*real(Slab(:,2));

%%% Plot an example of the z-axis.
z_x2a = real(Slab(ic2,2)) - 0.4*H;
z_y2a = m(ic2)*z_x2a + y_int(ic2);
p1 = [real(Slab(ic2,2)), imag(Slab(ic2,2))];
p2 = [z_x2a,  z_y2a];
dp_z = p2 - p1;

%%% Create the arrow
Arrow_Z = annotation('arrow', 'position', [p1(1), p1(2)+2*dp_z(2), dp_z(1), -dp_z(2)],...
    'linestyle','none');
set(Arrow_Z,'parent', gca);

%%% Create the line.
Line_Z = annotation('arrow', 'position', [p1(1), p1(2), dp_z(1), dp_z(2)]);
Line_Z.HeadStyle = 'none';
Line_Z.LineWidth = 2;
set(Line_Z,'parent', gca);

%%% Axis label.
text(p2(1)+0.5, p2(2), '$z$', 'FontSize', 20, 'Interpreter', 'latex')

%%%------------------------------------------------------------------------------------%%%
%%% Plot an example of the s-axis. Slope and y-intercept of tanget line.
m = tan(theta(ic2));
y_int = imag(Slab(ic2,2)) - m.*real(Slab(ic2,2));
z_y2b = imag(Slab(ic2,2)) + 0.4*H;
z_x2b = (1/m)*(z_y2b - y_int);
p3 = [z_x2b, z_y2b];
dp_s = p3 - p1;

%%% Create the arrow
Arrow_S = annotation('arrow', 'position', [p1(1), p1(2)+2*dp_s(2), dp_s(1), -dp_s(2)],...
    'linestyle','none');
set(Arrow_S,'parent', gca);

%%% Create the line.
Line_S = annotation('arrow', 'position', [p1(1), p1(2), dp_s(1), dp_s(2)]);
Line_S.HeadStyle = 'none';
Line_S.LineWidth = 2;
set(Line_S,'parent', gca);

%%% Axis label.
text(p3(1), p3(2)-0.5, '$s$', 'FontSize', 20, 'Interpreter', 'latex')

%%%------------------------------------------------------------------------------------%%%
%%% Add annotations for slab and shell thicknesses.
Arrow_Slab = annotation('doublearrow');
Arrow_Slab.Parent = gca;
Arrow_Slab.Position = [-5 0 0 H];
Arrow_Slab.LineWidth = 2;
Arrow_Slab.Head1Length = -10;
Arrow_Slab.Head2Length = -10;
text(-4.5, H/2, '$H$', 'FontSize', 20, 'Interpreter', 'latex')

Arrow_Shell = annotation('doublearrow');
Arrow_Shell.Parent = gca;
Arrow_Shell.Position = [-7.5 0 0 H_shell];
Arrow_Shell.LineWidth = 2;
Arrow_Shell.Head1Length = -10;
Arrow_Shell.Head2Length = -10;
text(-7, H_shell/2, '$H_{shell}$', 'FontSize', 20, 'Interpreter', 'latex')

%%%------------------------------------------------------------------------------------%%%
%%% Add annotiotions for the temperature conditions.
text(0.8*x2, H+1, '$T = T_b$', 'FontSize', 20, 'Interpreter', 'latex')
text(0.8*x2, -0.75, '$T = T_s$', 'FontSize', 20, 'Interpreter', 'latex')
chr = ['$T = T_b$' newline 'warm, convecting ice'];
text(0*x1, 0.7*H_shell, chr, 'FontSize', 18, 'Interpreter', 'latex')

%%%------------------------------------------------------------------------------------%%%
%%% Limits and such.
axis equal
ax = gca;
% XX = ax.XLim;
% ax.XLim = [-1 12];
% ax.YLim = [0 22];
ax.YDir = "reverse";
hold off
ax.YTickLabel = '';
ax.XTickLabel = '';
xlabel('Distance, $x$-axis', 'Interpreter', 'latex')
ylabel('Depth, $y$-axis', 'Interpreter', 'latex')


axesH = findall(fig, "Type", "axes");
set(axesH, "TickLabelInterpreter", 'latex')

%%%------------------------------------------------------------------------------------%%%
%%% Previous code
%%%------------------------------------------------------------------------------------%%%
%%% Plot the slab boundaries.
% plot([x1; real(Slab(:,1))], imag([Slab(1,1); Slab(:,1)]), 'k', LineWidth = 3)
% hold on
% plot([x1; real(Slab(:,2))], imag([Slab(1,2); Slab(:,2)]), 'k--', LineWidth = 2)
% plot([x1; real(Slab(:,3))], imag([Slab(1,3); Slab(:,3)]), 'k', LineWidth = 3)

%%% Plot the shell boundaries.
% plot([x1 x2], [0 0], 'k', LineWidth = 3)
% plot([x1 x2], x2*[1 1], 'k', LineWidth = 3)


%%% Find the depth where the plate interface hits the bottom of the conductive layer.
% [~, i_int] = min(abs(imag(Slab(:,3)) - H));
% plot([real(Slab(i_int,3)), x2], imag(Slab(i_int,3))*[1 1], 'k', LineWidth=3)

% Text_X = annotation('textbox') ;
% Text_X.Parent = gca;
% Text_X.Position = [10, 0, 2, 2];
% Text_X.FontSize = 20;
% Text_X.Interpreter = 'latex';
% Text_X.String = '$x$';
% Text_X.EdgeColor = 'none';

% Text_Y = annotation('textbox') ;
% Text_Y.Parent = gca;
% Text_Y.Position = [-0.5, 7.5, 2, 2];
% Text_Y.FontSize = 20;
% Text_Y.Interpreter = 'latex';
% Text_Y.String = '$y$';
% Text_Y.EdgeColor = 'none';

% %%% Plot an example of the s-axis and z-axis.
% z_x2a = real(Slab(ic2,2)) - H/4;
% z_y2a = m(ic2)*z_x2a + y_int(ic2);
% % plot([real(Slab(ic2,2)), z_x2a], [imag(Slab(ic2,2)), z_y2a], 'r', LineWidth=2)
% 
% p1 = [real(Slab(ic2,2)), imag(Slab(ic2,2))];
% p2 = [z_x2a,  z_y2a];
% dp = p2 - p1;
% quiver(p1(1), p1(2), dp(1), dp(2), 0, LineWidth=2, Color='r')
% % annotation('arrow', [real(Slab(ic2,2)), z_x2a], [imag(Slab(ic2,2)), z_y2a])
% 
% %%% Slope and y-intercept of tanget line.
% m = tan(theta(ic2));
% y_int = imag(Slab(ic2,2)) - m.*real(Slab(ic2,2));
% z_y2b = imag(Slab(ic2,2)) + H/4;
% z_x2b = (1/m)*(z_y2b - y_int);
% p2 = [z_x2b, z_y2b];
% dp = p2 - p1;
% quiver(p1(1), p1(2), dp(1), dp(2), 0, LineWidth=2, Color='r')
% % plot([real(Slab(ic2,2)), z_x2b], [imag(Slab(ic2,2)), z_y2b], 'r', LineWidth=2)

% Arrow_Z = annotation('arrow') ;
% Arrow_Z.Parent = gca;
% Arrow_Z.Position = [p1(1), p1(2), dp_z(1), dp_z(2)];
% Arrow_Z.LineWidth = 2;

% Arrow_S = annotation('arrow') ;
% Arrow_S.Parent = gca;
% Arrow_S.Position = [p1(1), p1(2), dp_s(1), dp_s(2)];
% Arrow_S.LineWidth = 2;