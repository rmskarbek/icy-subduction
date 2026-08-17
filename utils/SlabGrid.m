function [x_Slab, y_Slab] = SlabGrid(Out)

%%% This function generates a mesh grid of the slab in the Cartesian coordinates
%%% system that is used to define the entire model geometry. The mesh grid can be 
%%% used for making contour plots of simulation output that show the curved slab 
%%% geometry, instead of a rectangular geometry.

%%%-------------------------------------------------------------------------------%%%
%%% Input.
%%%-------------------------------------------------------------------------------%%%
% Out      - A structure of simulation output generated from IcySubduction.m.

%%%-------------------------------------------------------------------------------%%%
%%% OUTPUT.
%%%-------------------------------------------------------------------------------%%%
% x_Slab - An array of x-coordinate values for each location in the slab.

% y_Slab - An array of y-coordinate values for each location in the slab.
%%%-------------------------------------------------------------------------------%%%

%%% Get the number of grid points used for the simulation. This corresponds to the
%%% number of grid points along the z-axis.
nz_Points = Out.p.Numerical.GridPoints;

%%% Get the Cartesian coordinates of the top and bottom of the slab.
x_SlabBottom = real(Out.Slab.Bottom);
y_SlabBottom = imag(Out.Slab.Bottom);

x_SlabTop = real(Out.Slab.Top);
y_SlabTop = imag(Out.Slab.Top);

%%% Get the number of points along the s-axis.
ns_Points = numel(x_SlabTop);

%%% Allocate arrays to store the slab mesh grid.
x_Slab = nan(nz_Points, ns_Points);
y_Slab = nan(nz_Points, ns_Points);

%%% Loop through the number of s-axis points to populate the slab coordinate arrays.
for i = 1:ns_Points
    x_Slab(:, i) = linspace(x_SlabTop(i,1), x_SlabBottom(i,1), nz_Points)';
    y_Slab(:, i) = linspace(y_SlabTop(i,1), y_SlabBottom(i,1), nz_Points)';
end
