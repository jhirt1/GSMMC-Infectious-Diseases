function [S, V] = S_to_V(S, V, p)
% This function converts an susceptible individual to vaccinated after by
% checking the vaccination probability
% p: is the vaccination probability

% Find currently susceptible individuals
idxv = find(S == 1);
if ~isempty(idxv)
    %Roll a random number for each susceptible individual
    rolled_vac = rand(length(idxv), 1) < p;
    % update states
    V(rolled_vac) = 1;
    S(rolled_vac) = 0;
    %V(idxv) = rolled_vac;
    %S(idxv) = S(idxv) - rolled_vac;

end