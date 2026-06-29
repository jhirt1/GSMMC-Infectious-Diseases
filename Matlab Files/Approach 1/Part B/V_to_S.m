function [V,S] = V_to_S(V,S,VID,SAVdays)
% This funtion moves vaccinated individuals back to susceptible after
% vaccination immunity duration (VID)

% Find vaccinated individuals
ind = find(V);

% Check whose vaccine immunity has expired
vac = SAVdays(ind) >= VID;

% Move them back to susceptible
S(ind) = S(ind) | vac;

% Remove them from vaccinated class
V(ind) = V(ind)-vac;

end