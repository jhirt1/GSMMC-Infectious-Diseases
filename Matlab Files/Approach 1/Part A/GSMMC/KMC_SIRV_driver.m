% A simple KMC simulation of an SIR model
% June 08, 2026
tic
%clear; close all; clc
rng(10); % seed the random number generator
N = 100;  % size of the population
rho = 0.5;  % density of population per unit area
p = 0.01;

L = sqrt(N/rho); % length of square simulation area
sig2 = (0.025*L)^2; % infection varies over the whole length
rcd = 14; % It takes rcd days for an infected person to recover 
          % Recovery days can be a probability distribution
rsd = 1; % Number days of people stayed after recovery before becoming susceptible
vdays = 1; % Number days people after recovery before getting vaccinated
ni_val = [1, 10];   % initial number of infected individuals
sA_val1 = rand(N,1);
sA_val2 = sA_val1;
sA_val2(1:N/2) = 0 * sA_val2(1:N/2);
sA_val2 = sA_val2(randperm(N));
sA = sA_val2;
%why sA is defined as a vector of ones?
nv = 0; % initial number of vaccinated individuals
tspan = 365; % Run for 720 time units
% Initialize arrays to store results


nsim = 500; % average of 10 simulations (this needs to be much higher)
for init_ni = 1:length(ni_val)
    ni = ni_val(init_ni);
    % for sAact = 1:size(sA_val,2)
    %     sA = sA_val(:, sAact);
    %     ratios = sA ./ sA_val1;
    % 
    %     % 2. Find rows where the ratio is NOT 1 (using a tiny tolerance for decimals)
    %     changed_indices = find(abs(ratios - 1) > 1e-5);
    %     num_changed = length(changed_indices);
    %     prcnt = 100 * num_changed / N;

        susceptible = zeros(tspan, 1);
        infected = zeros(tspan, 1);
        recovered = zeros(tspan, 1);
        vaccinated = zeros(tspan,1);
        for sim = 1:nsim
            % Initialize health status for each simulation
            [S, H, R, V, Hdays,SARdays, sA, x, y] = HealthStatus_V(N, L, ni, sA); % Set individual health characteristics
            ImmunityTime = zeros(N,1);
            for t = 1:tspan
                % Update health status based on the SIR model dynamics
                %[S, H, R, V, x, y] = updateHealthStatus_V(S, H, R, V, Hdays, SARdays, sA, x, y,sig2,rcd,rsd, vdays);
                %[S, H, R, V, x, y, ImmunityTime] = updateHealthStatus_V(S, H, R, V, Hdays, SARdays, sA, x, y,sig2,rcd, vdays, ImmunityTime);
                [S, H, R, V, x, y, ImmunityTime] = updateHealthStatus_V(S, H, R, V, Hdays, SARdays, sA, x, y,sig2,rcd, ImmunityTime, p);
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
    figure;
    plot(ttime, susceptible,'b-', ttime, infected, 'r-.', ttime, recovered,'g:', ttime, vaccinated,'k-', LineWidth=2)
    legend('S','I','R', 'V',Location='northwest');
    % title(sprintf('Num of initial infected people = %d, with baseline activity', ni))
    % title(sprintf('Num of initial infected people = %d, 25%% of population has 0.5x activity', ni))
    title(sprintf('Num of initial infected people = %d, 50%% of population has 0 activity', ni))

    xlabel('hours')
    % end  
end   
toc







% %% Animate S-I-R-V states
% 
% Nanim = round(susceptible(1) + infected(1) + recovered(1) + vaccinated(1));
% 
% x = rand(Nanim,1);
% y = rand(Nanim,1);
% 
% figure
% 
% for k = 1:tspan
% 
%     clf
% 
%     % Current compartment counts
%     S = round(susceptible(k));
%     I = round(infected(k));
%     R = round(recovered(k));
% 
%     % Force total population to remain Nanim
%     V = Nanim - (S + I + R);
% 
%     % Safety check
%     S = max(S,0);
%     I = max(I,0);
%     R = max(R,0);
%     V = max(V,0);
% 
%     % Randomly assign individuals to states
%     idx = randperm(Nanim);
% 
%     state = zeros(Nanim,1);
% 
%     state(idx(1:S)) = 1;
%     state(idx(S+1:S+I)) = 2;
%     state(idx(S+I+1:S+I+R)) = 3;
%     state(idx(S+I+R+1:S+I+R+V)) = 4;
% 
%     hold on
% 
%     % Susceptible (blue)
%     scatter(x(state==1),y(state==1),100,[0.2 0.4 0.9],'filled')
% 
%     % Infected (red)
%     scatter(x(state==2),y(state==2),100,'r','filled')
% 
%     % Recovered (green)
%     scatter(x(state==3),y(state==3),100,[0 0.6 0],'filled')
% 
%     % Vaccinated (black)
%     scatter(x(state==4),y(state==4),100,'k','filled')
% 
%     axis([0 1 0 1])
%     axis square
%     axis off
% 
%     title(sprintf('Day %d',k))
% 
%     drawnow
%     pause(0.1)
% 
% end






v = VideoWriter('kmc_simulation.mp4','MPEG-4');
v.FrameRate = 10;
open(v);

for k = 1:tspan

    clf
    % ---- your plotting code here ----
    scatter(x(state==1),y(state==1),100,[0.2 0.4 0.9],'filled')
    scatter(x(state==2),y(state==2),100,'r','filled')
    scatter(x(state==3),y(state==3),100,[0 0.6 0],'filled')
    scatter(x(state==4),y(state==4),100,'k','filled')

    axis([0 1 0 1])
    axis square
    axis off

    title(sprintf('Day %d',k))

    drawnow

    % capture frame
    frame = getframe(gcf);
    writeVideo(v, frame);

end

close(v);