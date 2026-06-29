function [Hdays] = updateInfectionDays(H,Hdays)

ind=find(H); % Find the infected individuals
Hdays(ind) = Hdays(ind) + 1;