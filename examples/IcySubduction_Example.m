%%% This script will run a Europa subduction simulation; produce some plots of the
%%% simulation output; analyze the output to compute the bending, buoyancy, and shear
%%% resistance forces; and plot the forces.

%%% Generate a parameters structure with all necessary information to run a 
%%% simulation.
p = IcySubduction_Parameters;

%%% Run the simulation.
Out = IcySubduction(p);

%%% Plots of the temperature, porosity, and bulk density fields.
IcySubductionPlots(Out);

%%% Calculate the forces associated with plate bending.
Bend = BendingForce(Out);

%%% Calculate the force due to shear resistance along the plate interface.
% mu_flag = 'constant';
mu_flag = 'temperature';
Shear = IceShearResistance(Out, mu_flag);

%%% Calculate the buoyancy force.
Buoyancy = BuoyancyForce(Out);

%%% Plot some failure envelopes.
FailureEnvelopesPlot(Bend, Out)

%%% Plot the bending, shear, and buoyancy forces.
ForcesPlot(Bend, Buoyancy, Shear, Out);