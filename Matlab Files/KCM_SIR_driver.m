% A simple KMC simulation of an SIR model
% June 08, 2026
clear;
rng(10); % seed the random number generator
N = 50;  % size of the population
rho = 0.5;  % density of population per unit area

L = sqrt(N/rho); % length of square simulation area
sig2 = (0.025*L)^2; % infection varies over the whole length
rcd = 50; % It takes rcd days for an infected person to recover 
          % Recovery days can be a probability distribution
ni = 1;   % initial number of infected individuals

tspan = 720; % Run for 720 time units
% Initialize arrays to store results
susceptible = zeros(tspan, 1);
infected = zeros(tspan, 1);
recovered = zeros(tspan, 1);

nsim = 500; % average of 10 simulations (this needs to be much higher)
for sim = 1:nsim
    % Initialize health status for each simulation
    [S, H, R, Hdays, sA, x, y] = HealthStatus(N, L, ni); % Set individual health characteristics
    
    for t = 1:tspan
        % Update health status based on the SIR model dynamics
        [S, H, R, x, y] = updateHealthStatus(S, H, R, Hdays, sA, x, y,sig2,rcd);
        % Update number of days of infection
        Hdays = updateInfectionDays(H,Hdays);
        [x,y] = updatePositions(x,y,L,N,sig2);

        % Store results
        susceptible(t) = susceptible(t) + sum(S);
        infected(t) = infected(t) + sum(H);
        recovered(t) = recovered(t) + sum(R);
    end
end

ttime = 1:tspan;
susceptible = susceptible/nsim;
infected = infected/nsim;
recovered = recovered/nsim;


plot(ttime, susceptible,'b-', ttime, infected, 'r-.', ttime, recovered,'g:',LineWidth=2)
legend('S','I','R',Location='northwest');
xlabel('hours')