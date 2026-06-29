function [S, H] = SusceptibleToInfected(S,H,sA,x,y,sig2)
% Function to determine if susceptible individuals become infected

tsig2=2*sig2;
% Find indices of susceptible individuals
indi = find(S); 
% Find indices of infected individuals
indj = find(H); 

% 
for id=1:length(indi) % Loop over infected individuals
    i=indi(id); 
    xi=x(i); yi=y(i); pi=0;
    for jd=1:length(indj)
        j=indj(jd);
        rij2=(xi-x(j))^2 + (yi-y(j))^2;
        pi = pi + sA(j)*exp(-rij2/tsig2);
    end
    ii = floor(pi + rand);
    if ii
        H(i) = 1;
        S(i) = 0;
    end
end