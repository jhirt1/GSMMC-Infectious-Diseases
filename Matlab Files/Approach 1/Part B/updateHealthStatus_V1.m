function [S, H, R, V, x, y] = updateHealthStatus_V1(S, H, R, V, Hdays, SARdays, SAVdays, sA, x, y,sig2,rcd, ID,VID,p)
%function [S, H, R, V,ImmunityTime, x, y] = updateHealthStatus_V1(S, H, R, V, Hdays, SARdays, SAVdays, sA, x, y,sig2,rcd, ImmunityTime,VID,p)

[S, H] = SusceptibleToInfected(S,H,sA,x,y,sig2);
%[H, R, ImmunityTime] = InfectedToRecovered_V(H, R, Hdays, rcd, ImmunityTime);
%[R,S,SARdays] = RecoveredToSusceptible_V(R,S,SARdays,ImmunityTime);
[H, R, SARdays] = InfectedToRecovered(H,R,Hdays,rcd, SARdays);
[R, S] = R_to_S(R,S,ID,SARdays);
[S,V] = S_to_V(S,V,p);
[V,S] = V_to_S(V,S,VID,SAVdays);

