function dudt = fSIR(t,u,beta,gamma)
    S = u(1); I=u(2); R=u(3);
    dudt = [-beta*S*I; beta*S*I - gamma*I; gamma*I];
end