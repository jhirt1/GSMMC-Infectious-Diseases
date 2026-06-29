%My original function 
function [S, V] = SusceptibleToVaccinated(S, V, p)
    % Find currently susceptible individuals
    idxv = find(S == 1); 

    if ~isempty(idxv)
        % Roll a random number for each susceptible individual
        rolled_vac =  rand(length(idxv), 1) < p; 
        %length(rolled_vac)
        %V(rolled_vac)
        %pause
        %idx_to_vaccinate = idxv(rolled_vac);

        % Update states
        % V(rolled_vac) = 1;
        % S(rolled_vac) = 0;
        V(idxv) = rolled_vac;
        S(idxv) = S(idxv) - rolled_vac;
    end
end