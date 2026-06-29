function [R, S] = R_to_S(R,S,ID,SARdays)
% this function moves people from recovered to susceptible
ind = find(R); % Find the susceptible individuals
% Check who has stayed recovered long enough
rec = SARdays(ind) >= ID;
    
% Add to S (don't overwrite — use OR logic)
S(ind) = S(ind) | rec;

% Remove from R
R(ind) = R(ind) - rec;