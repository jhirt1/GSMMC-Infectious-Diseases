function [S,H,R,Hdays,sA,x,y] = HealthStatus(N,L,ni)
%===========================================
% Function to assign initial positions and health status
% Input
% N: population size
% L: length of square area
% ni: size of initial infected population
% Output
% S: indicator for suceptible individuals
% H: indicator for infected individuals
% R: indicator for recovered individuals
%sA: Social activity of a person
% (x,y) : x,y coordinates of individuals
%===========================================

% Initialize the positions of individuals randomly within the area
x = rand(N, 1) * L; % x coordinates 
y = rand(N, 1) * L; % corresponding y coordinates

% Preallocate arrays for each characteristic of an individual
S = zeros(N,1); % susceptible individuals: 1 - susceptible, 0 - not susceptible
% Initialize all individuals as susceptible
S(:) = 1;  % Set all individuals to susceptible (1)

H = zeros(N,1); % infected population: 1 - infected, 0 - not infected
% number of days an individual has been infected.  
Hdays = zeros(N,1); % Initialized to zero
% Pick randomly in the range [1,N]
ind = randperm(N,ni); % pick ni individuals randomly from N to be infected

H(ind) = 1; % Set the selected individuals to infected
S(ind) = 0; % So the individual is not susceptible

R = zeros(N,1); % initial recovered population: 1 -recovered, 0 - not

sA=ones(N,1); %rand(N,1);