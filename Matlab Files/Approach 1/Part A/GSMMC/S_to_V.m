function [S, V] = S_to_V(S,V, H, nu0, coverage)

% This function converts an susceptible individual to vaccinated after by
% checking the vaccination rate and the vaccination efficacy
% Population size
N = length(S);
% Current vaccination coverage
currentCoverage = sum(V)/N;

% Vaccination stops if maximum coverage has been reached
if currentCoverage < coverage

    % Infection prevalence
    prevalence = sum(H)/N;

    % Vaccination probability
    nu = nu0*prevalence^2;

    % Determine who gets vaccinated
    Vnew = S & (rand(N,1) < nu);
    S = S - Vnew;
    V = V + Vnew;

end