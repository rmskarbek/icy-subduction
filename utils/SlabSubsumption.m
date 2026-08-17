function [Out, ArcLengthSubsumed, DepthSubsumed, i_Subsumed, T_subsumed]...
    = SlabSubsumption(Out)

%%% This function finds locations in the slab where the temperature is equal to the
%%% temperature in the convecting ice. This serves as a proxy for material that has
%%% been subsumed. The slab is fully subsumed when at the arc length where the
%%% temperature of the entire column is greater than or equal to the convecting ice
%%% temperature.

%%%-------------------------------------------------------------------------------%%%
%%% INPUT.
%%%-------------------------------------------------------------------------------%%%
% Out      - A structure of simulation output generated from IcySubduction.m.


%%%-------------------------------------------------------------------------------%%%
%%% OUTPUT.
%%%-------------------------------------------------------------------------------%%%
% Out               - The arc length and depth of subsumption are added to Out.Slab.

% ArcLengthSubsumed - The arc length along the center line of the slab where
%                     subsumption occurs.

% DepthSubsumed     - The depth at the bottom of the column when subsumption occurs,
%                     minus the slab thickness.

% i_subsumed        - The array column index corresponding to ArcLengthSubsumed.

% T_subsumed        - The temperature used to determine where subsumption occurs.
%%%-------------------------------------------------------------------------------%%%

%%% Get the temperature simulation data and define the subsumption temperature.
T_subsumed = Out.p.Temperature.BasalTemp - 0.1;
% T_subsumed = Out.p.Temperature.BasalTemp - 1;

Temp_min = min(Out.Temperature,[],1)';
i_Subsumed = find(Temp_min > T_subsumed, 1);

%%% If the slab was not fully subsumed, set the index and arc length equal to the
%%% full slab length.
if isempty(i_Subsumed)
    i_Subsumed = numel(Out.Slab.ArcLength);
    ArcLengthSubsumed = nan;
    DepthSubsumed = nan;
else
    ArcLengthSubsumed = Out.Slab.ArcLength(i_Subsumed);
    DepthSubsumed = Out.Slab.Depth(end, i_Subsumed) - Out.p.Geometry.SlabThick;
end

Out.Slab.ArcLengthSubsumed = ArcLengthSubsumed;
Out.Slab.DepthSubsumed = DepthSubsumed;