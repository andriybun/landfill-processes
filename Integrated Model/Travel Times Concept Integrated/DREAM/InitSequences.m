function [Sequences] = InitSequences(X,Sequences,MCMCPar);
% Initialize sequences

for qq  = 1:MCMCPar.seq,
    % Initialize Sequences
    Sequences(1,1:MCMCPar.n+2,qq) = X(qq,1:MCMCPar.n+2);
end