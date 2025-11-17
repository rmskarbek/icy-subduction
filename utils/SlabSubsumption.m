function [ArcLengthSubsumed, i_Subsumed] = SlabSubsumption(Out)

%%% This function finds locations in the slab where the temperature is equal to the
%%% temperature in the convecting ice. This serves as a proxy for material that has
%%% been subsumed. The slab is fully subsumed when at the arc length where the
%%% temperature of the entire column is greater than or equal to the convecting ice
%%% temperature.

%%% Get the temperature simulation data and define the subsumption temperature.
T_subsumed = Out.p.Temperature.BasalTemp - 0.01;

Temp_min = min(Out.Temperature,[],1)';
i_Subsumed = find(Temp_min > T_subsumed, 1);

%%% If the slab was not fully subsumed, set the index and arc length equal to the
%%% full slab length.
if isempty(i_Subsumed)
    i_Subsumed = numel(Out.Slab.ArcLength);
end
ArcLengthSubsumed = Out.Slab.ArcLength(i_Subsumed);


%%%-------------------------------------------------------------------------------%%%
%%% This block of code finds the arc length where the slab begins to subsume, rather
%%% than where it is fully subsumed.

% %%% Find locations where the slab temperature is equal to or greater than the
% %%% subsumption temperature, and replace those temperature values with nan.
% Temp = Out.Temperature;
% i_subsumed = Temp >= T_subsumed;
% Temp(i_subsumed) = nan;
% 
% %%% Check if the slab was fully subsumed, and if so determine the arc length where
% %%% that occurs. First remove the bottom row of temperature values, because that is 
% %%% where the basal boundary condition is set.
% Temp = Temp(1:end-1,:);
% ArcLengthSubsumed = nan;
% 
% i_Subsumed = 1;
% while sum(isnan(Temp(:,i_Subsumed))) == 0
%     i_Subsumed = i_Subsumed + 1;
% end
% 
% if i_Subsumed < size(Temp,2)
%     ArcLengthSubsumed = Out.Slab.ArcLength(i_Subsumed);
% end