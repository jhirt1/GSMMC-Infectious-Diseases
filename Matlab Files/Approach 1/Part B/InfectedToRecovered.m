function [H, R, SARdays] = InfectedToRecovered(H,R,Hdays,rcd, SARdays)
% This function converts an infected individual to recovered after
% rcd days of infection. The individual does not join the susceptible pool

ind = find(H); % Find the infected individuals
rec=Hdays(ind) == rcd; % If infection days = rcd...
R(ind) = rec; %...set to recovered
H(ind) = H(ind)-rec; % set infection to zero
SARdays(ind) = 1; 
