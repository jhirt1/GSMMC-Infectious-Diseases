function dudt = fSIRV(t,u,beta,gamma,nu,p)
    S = u(1); I=u(2); R=u(3); V=u(4);
    betaSI = beta*S*I; nuR = nu*R; pS = p*S; gammaI=gamma*I;
    dudt = [-betaSI+nuR-pS; betaSI - gammaI; gammaI-nuR; pS];
end