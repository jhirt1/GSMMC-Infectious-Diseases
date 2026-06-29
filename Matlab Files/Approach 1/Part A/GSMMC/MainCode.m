% A simple KMC simulation of an SIR model
% June 08, 2026
clear all; clc
rng(10); % seed the random number generator
N = 50;  % size of the population
rho = 0.5;  % density of population per unit area

L = sqrt(N/rho); % length of square simulation area
sig2 = (0.025*L)^2; % infection varies over the whole length
rcd = 50; % It takes rcd days for an infected person to recover 
          % Recovery days can be a probability distribution
rsd = 1; % Number days of people stayed after recovery before becoming susceptible
%vdays = round(normrnd(5, 1));
% vdays = max(0, vdays);
vdays=1;
% Number days people after recovery before getting vaccinated
ni = 1;   % initial number of infected individuals
nv = 0; % initial number of vaccinated individuals
tspan = 720; % Run for 720 time units
% vdays = round(normrnd(20, 10, [N, 1]));
%vdays = 1;

% Initialize arrays to store results
susceptible = zeros(tspan, 1);
infected = zeros(tspan, 1);
recovered = zeros(tspan, 1);
vaccinated = zeros(tspan,1);

nsim = 50; % average of 10 simulations (this needs to be much higher)
for sim = 1:nsim
    % Initialize health status for each simulation
    [S, H, R, V, Hdays,SARdays, sA, x, y] = HealthStatus_V(N, L, ni); % Set individual health characteristics
    
    for t = 1:tspan
        % Update health status based on the SIR model dynamics
        [S, H, R, V, x, y] = updateHealthStatus_V(S, H, R, V, Hdays, SARdays, sA, x, y,sig2,rcd,rsd, vdays);
        % Update number of days of infection and recovery
        Hdays = updateInfectionDays(H,Hdays);
        SARdays = updateRecoveryDays(R,SARdays);
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


%plot(ttime, susceptible,'b-', ttime, infected, 'r-.', ttime, recovered,'g:',LineWidth=2)
%legend('S','I','R',Location='northwest');
%xlabel('hours')


plot(ttime, susceptible,'b-', ttime, infected, 'r-.', ttime, recovered,'g:', ttime, vaccinated,'k-', LineWidth=2)
legend('S','I','R', 'V',Location='northwest');
xlabel('hours')