% A simple KMC simulation of an SIR model
% June 08, 2026
clear; close all; clc;
rng(10); % seed the random number generator
N = 1000;  % size of the population
rho = 0.5;  % density of population per unit area

L = sqrt(N/rho); % length of square simulation area
sig2 = (0.025*L)^2; % infection varies over the whole length
rcd = 14; % It takes rcd days for an infected person to recover 
          % Recovery days can be a probability distribution

p = 0.008; % vaccination probability for python 
%p = 0.01; 
ni = 10;   % initial number of infected individuals

ID = 45; % natural  immunity duration
VID = 120; % Vaccine immunity duration

tspan = 365; % Run for tspan time units

% Initialize arrays to store results
susceptible = zeros(tspan, 1);
infected = zeros(tspan, 1);
recovered = zeros(tspan, 1);
vaccinated = zeros(tspan, 1);

nsim = 500; % average of nsim simulations (this needs to be much higher)
for sim = 1:nsim
    % Initialize health status for each simulation
    [S, H, R,V, Hdays,SARdays, SAVdays,  sA, x, y] = HealthStatus_V(N, L, ni); % Set individual health characteristics
    %ImmunityTime = zeros(N,1);
    for t = 1:tspan
        

        % Update health status based on the SIR model dynamics
        [S, H, R, V, x, y] = updateHealthStatus_V1(S, H, R, V, Hdays, SARdays, SAVdays, sA, x, y,sig2,rcd, ID,VID,p);
        %[S, H, R, V,ImmunityTime, x, y] = updateHealthStatus_V1(S, H, R, V, Hdays, SARdays, SAVdays, sA, x, y,sig2,rcd, ImmunityTime,VID,p);
        % Update number of days of infection, days of recovery, and days of
        % vaccinated
        Hdays = updateInfectionDays(H,Hdays);
        SARdays = updateRecoveryDays(R,SARdays);
        SAVdays = updateVaccinatedDays(V, SAVdays);
        [x,y] = updatePositions(x,y,L,N,sig2);

        % Store results
        susceptible(t) = susceptible(t) + sum(S);
        infected(t) = infected(t) + sum(H);
        recovered(t) = recovered(t) + sum(R);
        vaccinated(t) = vaccinated(t) + sum(V);
    end
end

ttime = 1:tspan;
susceptible = susceptible/nsim;
infected = infected/nsim;
recovered = recovered/nsim;
vaccinated = vaccinated/nsim;


plot(ttime, susceptible,'b-', ttime, infected, 'r-.', ttime, recovered,'g:', ttime, vaccinated,'k-', LineWidth=2)
legend('S','I','R', 'V',Location='northwest');
xlabel('days')
ylabel('population size')
%title(sprintf('SIRVS KMC  |  N=%d, T=%d, I_0=%d, VID=%d, H_{days}=%d', ...
    %N, tspan, ni, VID, rcd));
title(sprintf('SIRVS KMC  |  N=%d, T=%d, I_0=%d, ID=%d, VID=%d, H_{days}=%d', ...
   N, tspan, ni, ID, VID, rcd));