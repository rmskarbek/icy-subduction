%%% This script will run a Europa subduction simulation; produce some plots of the
%%% simulation output; analyze the output to compute the bending and frictional forces;
%%% and plot the forces.

%%% Generate a parameters structure will all necessary information to run a simulation.
p = IcySubduction_Parameters;

%%% Schematic plot of the model geometry with the slab, ice shell and coordinate axes.
%%% NOTE: This code is optimized for H = 5 km and H_shell = 25 km, but should still work
%%% for other thicknesses.
GeometryPlot(p);

%%% Run the simulation.
Out = IcySubduction(p);

%%% Plots of the temperature, viscosity, and porosity fields.
IcySubduction_Plots1(Out);

%%% Calculate the forces associated with plate bending.
Bend = BendingForce(Out);

%%% Calculate the force due to shear resistance along the plate interface.
% mu_flag = 'constant';
mu_flag = 'temperature';
Shear = IceShearResistance(Out, mu_flag);

%%% Plot some failure envelopes and plot the bending and shear forces.
ForcesPlot1(Bend, Shear, Out);