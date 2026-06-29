function [SAVdays] = updateVaccinatedDays(V,SAVdays)

ind=find(V); % Find the vaccinated individuals
SAVdays(ind) = SAVdays(ind) + 1;
