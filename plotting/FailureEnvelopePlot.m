%%% Run a subduction simulation.
p = IcySubduction_Parameters;
Out = IcySubduction(p);

%%% Compute the failure envelope.
Temp = Out.Temperature;
Sigma_L = 1e-6*Out.VerticalStress;
Depth = Out.Slab.Depth;
mu_flag = 'temperature';

[S_Compression, Tau_Compression, S_Tension, Tau_Tension, Sigma_Diff]...
    = IceFailureEnvelope(Depth(:,1), Temp(:,1), Sigma_L(:,1), mu_flag);

%%% Convert to km and MPa for plotting.
Depth = 1e-3*Depth;
S_Compression = 1e-6*S_Compression;
Tau_Compression = 1e-6*Tau_Compression;
S_Tension = 1e-6*S_Tension;
Tau_Tension = 1e-6*Tau_Tension;
Sigma_Diff = 1e-6*Sigma_Diff;

%%% Set up the figure.
font = 'Palatino Linotype';
figure('DefaultTextFontName',font,'DefaultAxesFontName',font,...
    'DefaultAxesFontSize', 16);

plot(Tau_Compression, Depth(:,1), 'g', LineWidth = 2)
hold on
plot(Sigma_Diff, Depth(:,1), 'm', LineWidth = 2)
plot(S_Compression, Depth(:,1), 'w--', LineWidth = 2)
hold off
ax = gca;
ax.YDir = "reverse";
ax.XLim = [0 max(Tau_Compression) + 1];
xlabel('Failure Stress (MPa)')
ylabel('Depth (km)')

legend(['Brittle Failure' newline '(friction)'], ['Ductile Failure' newline...
    '(flow laws)'], 'Failure Envelope')