function [x,y] = updatePositions(x,y,L,N,sig2)

sig=sqrt(sig2);
% Compute displacements as a Brownian motion
dx = normrnd(0, sig, [N, 1]); % normal distribution, mean 0, std = sig;
dy = normrnd(0, sig, [N, 1]);

% Update positions based on the computed displacements
x = mod(x + dx, L); % periodic boundary conditions 
y = mod(y + dy, L);
%i dont understand how x and y are changing
