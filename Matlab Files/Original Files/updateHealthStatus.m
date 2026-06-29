function [S, H, R, x, y] = updateHealthStatus(S, H, R, Hdays, sA, x, y,sig2,rcd)

[H, R] = InfectedToRecovered(H,R,Hdays,rcd);
[S, H] = SusceptibleToInfected(S,H,sA,x,y,sig2); 