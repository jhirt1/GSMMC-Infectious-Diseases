function [S,H,R,V,Hdays,SARdays,SAVdays, sA,x,y] = HealthStatus_V(N,L,ni)
%===========================================
% Function to assign initial positions and health status
% Input
% N: population size
% L: length of square area
% ni: size of initial infected population
% nv : size of initial vaccinated population
% Output
% S: indicator for suceptible individuals
% H: indicator for infected individuals
% R: indicator for recovered individuals
% V: indicator for vacinated individuals
% sA: Social activity of a person
% (x,y) : x,y coordinates of individuals
% Hdays: Number of days people have stayed infected
% SARdays: Number of days people have stayed recovered
% SAVdays: Number of days people have stayed vaccinated
%===========================================

% Initialize the positions of individuals randomly within the area
x = rand(N, 1) * L; % x coordinates 
y = rand(N, 1) * L; % corresponding y coordinates

% Preallocate arrays for each characteristic of an individual
S = zeros(N,1); % susceptible individuals: 1 - susceptible, 0 - not susceptible

% Initialize all individuals as susceptible
S(:) = 1;  % Set all individuals to susceptible (1)

H = zeros(N,1); % infected population: 1 - infected, 0 - not infected
% number of days an individual has been infected,recovered and vaccinated.  
Hdays = zeros(N,1); % Initialized to zero
SARdays = zeros(N,1); % Initialized to zero
SAVdays = zeros(N,1); % Initialized to zero

% Pick randomly in the range [1,N]
ind = randperm(N,ni); % pick ni individuals randomly from N to be infected

H(ind) = 1; % Set the selected individuals to infected
S(ind) = 0; % So the individual is not susceptible

R = zeros(N,1); % initial recovered population: 1 -recovered, 0 - not

V = zeros(N,1); % initial vaccinated population: 1 - vaccinated, 0 - not

sA=rand(N,1); %ones(N,1);