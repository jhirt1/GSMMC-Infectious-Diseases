function [S, H, R, V, x, y, ImmunityTime] = updateHealthStatus_V(S, H, R, V, Hdays, SARdays, sA, x, y,sig2,rcd, ImmunityTime, p)

% [H, R] = InfectedToRecovered(H,R,Hdays,rcd);
[H, R, ImmunityTime] = InfectedToRecovered(H,R,Hdays,rcd,ImmunityTime);
[S, H] = SusceptibleToInfected(S,H,sA,x,y,sig2); 
%[S, V] = SusceptibleToVaccinated(S,V, vdays, SARdays);
%[R, S] = RecoveredToSusceptible(R,S,rsd, SARdays);
%rsd : how many days you stay after recovery before you become
% susceptible.
% vdays the nuber of days people stay before the get vaccinated

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[R, S] = RecoveredToSusceptible(R,S,SARdays,ImmunityTime);
[S, V] = SusceptibleToVaccinated(S, V, p);