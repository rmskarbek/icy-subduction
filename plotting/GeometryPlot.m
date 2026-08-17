function GeometryPlot

%%%-------------------------------------------------------------------------------%%%
%%% This function plots the model geometry with the slab, ice shell and coordinate 
%%% axes. It also annotates the figure to the show the temperature boundary 
%%% conditions.

%%% NOTE: This code is optimized for H = 5 km and H_shell = 25 km, but should still 
%%% work for other thicknesses.
%%%-------------------------------------------------------------------------------%%%

%%% Get the geometry of the slab.
% p = IcySubduction_Parameters;
% H = p.Geometry.SlabThick;
% H_shell = p.Geometry.ShellThick;
% R_min = p.Geometry.CurveRadius;

H = 5e3;
H_shell = 25e3;
R_min = 2*H;

GeoFlag = 'Buffett';                % a flag that determines the slab geometry.

switch GeoFlag
    case 'Buffett'
        type = 'geometry';
        [w_top, w_center, w_bottom, theta, ~, ~, ~] = Buffett2006(type, H, H_shell,...
            R_min);

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
theme(fig,"light")
set(gcf,'Color','w')

%%% Limits for the x-axis.
xf = ceil(max(H_shell, max(real(Slab(:,3)))));
x1 = -10;
x2 = xf;

%%% Plot the convective layer.
area([x1 x2], H_shell*[1 1], FaceColor = [0.18 0.75 0.94]);
hold on

%%% Plot the conductive layer.
area([x1 x2], H*[1 1], FaceColor = '#0072BD');

%%% Plot the slab.
RGB = hex2rgb('#0072BD');
fill(real([x1; Slab(:,3); flip(Slab(:,1)); x1]), imag([0; Slab(:,3); flip(Slab(:,1));...
    1i*H]), RGB)

%%%-------------------------------------------------------------------------------%%%
%%% Plot cartesian coordinate axes.
Arrow_X = annotation('arrow') ;
Arrow_X.Parent = gca;
Arrow_X.Position = [0, 0, 10, 0];
Arrow_X.LineWidth = 2;
Arrow_X.Color = 'k';

Arrow_Y = annotation('arrow') ;
Arrow_Y.Parent = gca;
Arrow_Y.Position = [0, 0, 0, 10];
Arrow_Y.LineWidth = 2;
Arrow_Y.HeadLength = -10;
Arrow_Y.Color = 'k';

text(10, -0.5, '$x$', 'FontSize', 20, 'Interpreter', 'latex', 'Color', 'k')
text(0.5, 10, '$y$', 'FontSize', 20, 'Interpreter', 'latex', 'Color', 'k')

%%%-------------------------------------------------------------------------------%%%
%%% Plot the center line and the trench location.
plot(real(Slab(:,2)), imag(Slab(:,2)), 'k--', LineWidth = 1)

%%% Plot the plate interface.
[~, i_int] = min(abs(imag(Slab(:,3)) - H));
i_int = i_int + 1;
plot(real(Slab(1:i_int,3)), imag(Slab(1:i_int,3)), 'r--', LineWidth = 2)

% %%% Plot some locations of the column.
% plot(real(Slab(1, [1 3])), imag(Slab(1, [1 3])), 'k', LineWidth=3)
ic2 = round(0.3*numel(theta));
% ic2 = round(0.35*numel(theta));

% plot(real(Slab(ic2, [1 3])), imag(Slab(ic2, [1 3])), 'k', LineWidth=3)
ic3 = round(0.25*numel(theta));
plot(real(Slab(ic3, [1 3])), imag(Slab(ic3, [1 3])), 'r', LineWidth=2)


%%%-------------------------------------------------------------------------------%%%
%%% Plot an example of the z-axis. Tangent angle (psi) along the center line.
psi = pi/2 - theta;
%%% Slope of normal lines.
m = -tan(psi);
%%% y-intercept of normal lines.
y_int = imag(Slab(:,2)) - m.*real(Slab(:,2));

z_x2a = real(Slab(ic2,2)) - 0.4*H;
z_y2a = m(ic2)*z_x2a + y_int(ic2);
p1 = [real(Slab(ic2,2)), imag(Slab(ic2,2))];
p2 = [z_x2a,  z_y2a];
dp_z = p2 - p1;

%%% Create the arrow
Arrow_Z = annotation('arrow', 'position', [p1(1), p1(2)+2*dp_z(2), dp_z(1), -dp_z(2)],...
    'linestyle','none');
Arrow_Z.Color = 'k';
set(Arrow_Z,'parent', gca);

%%% Create the line.
Line_Z = annotation('arrow', 'position', [p1(1), p1(2), dp_z(1), dp_z(2)]);
Line_Z.HeadStyle = 'none';
Line_Z.LineWidth = 2;
Line_Z.Color = 'k';
set(Line_Z,'parent', gca);

%%% Axis label.
text(p2(1)+0.5, p2(2), '$z$', 'FontSize', 20, 'Interpreter', 'latex', 'Color', 'k')

%%%-------------------------------------------------------------------------------%%%
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
Arrow_S.Color = 'k';
set(Arrow_S,'parent', gca);

%%% Create the line.
Line_S = annotation('arrow', 'position', [p1(1), p1(2), dp_s(1), dp_s(2)]);
Line_S.HeadStyle = 'none';
Line_S.LineWidth = 2;
Line_S.Color = 'k';
set(Line_S,'parent', gca);

%%% Axis label.
text(p3(1), p3(2)-0.5, '$s$', 'FontSize', 20, 'Interpreter', 'latex', 'Color', 'k')

%%%-------------------------------------------------------------------------------%%%
%%% Plot the slab dip angle.
plot([p1(1) p1(1) + 0.75*H], p1(2)*[1 1], 'k', LineWidth=1)

r_theta = H/2;
theta_theta = linspace(theta(ic2), 0, 20)';
z_theta = r_theta*exp(1i*theta_theta);
plot(p1(1) + real(z_theta), p1(2) + imag(z_theta), 'k', LineWidth=1);
text(p1(1) + 0.5*H, p1(2) + 0.15*H, '$\theta(s)$', 'FontSize', 15,...
    'Interpreter', 'latex')


%%%-------------------------------------------------------------------------------%%%
%%% Add annotations for slab and shell thicknesses.
Arrow_Slab = annotation('doublearrow');
Arrow_Slab.Parent = gca;
Arrow_Slab.Position = [-8.25 0 0 H];
Arrow_Slab.LineWidth = 2;
Arrow_Slab.Head1Length = -10;
Arrow_Slab.Head2Length = -10;
Arrow_Slab.Color = 'k';
text(-9.75, H/2, '$H$', 'FontSize', 18, 'Interpreter', 'latex', 'Color', 'k')
% chr2 = ['$\rho_g = \rho_{slab}$' newline '$\phi = \phi(s,z)$'];
chr2 = ['$T = T(s,z)$' newline '$\phi = \phi(s,z)$'];
text(-7, H/2, chr2, 'FontSize', 18, 'Interpreter', 'latex') 

Arrow_Shell = annotation('doublearrow');
Arrow_Shell.Parent = gca;
Arrow_Shell.Position = [-7.5 0 0 H_shell];
Arrow_Shell.LineWidth = 2;
Arrow_Shell.Head1Length = -10;
Arrow_Shell.Head2Length = -10;
Arrow_Shell.Color = 'k';
text(-7, H_shell/2, '$H_{shell}$', 'FontSize', 18, 'Interpreter', 'latex', 'Color', 'k')

%%%-------------------------------------------------------------------------------%%%
%%% Add annotiotions for the temperature and density conditions.
text(0.75*x2, H+1, '$T = T_b$', 'FontSize', 18, 'Interpreter', 'latex')
text(0.75*x2, -0.75, '$T = T_s$', 'FontSize', 18, 'Interpreter', 'latex')
% chr3 = ['$\rho_g = \rho_{cond}$' newline '$\phi = \phi_0(y)$'];
chr3 = '$\phi = \phi_0(y)$';
text(0.75*x2, H/2, chr3, 'FontSize', 18, 'Interpreter', 'latex') 

% chr1 = ['$\rho_g = \rho_{conv}$' newline '$\phi = 0$' newline '$T = T_b$' newline...
%     'warm, convecting ice'];
chr1 = ['$T = T_b$' newline '$\phi = 0$' newline 'warm, convecting ice'];
text(0*x1, 0.7*H_shell, chr1, 'FontSize', 18, 'Interpreter', 'latex')

%%%-------------------------------------------------------------------------------%%%
%%% Limits and such.
axis equal
ax = gca;
ax.XLim = [x1 x2];
ax.YLim = [-2 H_shell];
ax.YDir = "reverse";
hold off
ax.YTickLabel = '';
ax.XTickLabel = '';
xlabel('Distance, $x$-axis', 'Interpreter', 'latex')
ylabel('Depth, $y$-axis', 'Interpreter', 'latex')


axesH = findall(fig, "Type", "axes");
set(axesH, "TickLabelInterpreter", 'latex')