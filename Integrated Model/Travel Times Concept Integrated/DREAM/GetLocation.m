function [x_old,p_old,log_p_old] = GetLocation(X,MCMCPar);
% Extracts the current location and density of the chain

% First get the current location
x_old = X(1:MCMCPar.seq,1:MCMCPar.n);

% Then get the current density
p_old = X(1:MCMCPar.seq,MCMCPar.n+1);

% Then get the current logdensity
log_p_old = X(1:MCMCPar.seq,MCMCPar.n+2);