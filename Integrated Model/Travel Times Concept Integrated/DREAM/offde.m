function [xnew,CR,alpha_s] = offde(X,Zoff,CR,MCMCPar,Update,Table_JumpRate,ParRange,BoundHandling)
% Calculate candidate points using discrete proposal distribution

% Determine how many pairs to use for each jump in each chain
[DEversion] = DEStrategy(MCMCPar);

% Generate uniform random numbers for each chain to determine which dimension to update
D = rand(MCMCPar.seq,MCMCPar.n);

% Generate noise to ensure ergodicity for each individual chain
noise_x = MCMCPar.eps * (2 * rand(MCMCPar.seq,MCMCPar.n) - 1);

% Initialize the delta update to zero
xjump = zeros(MCMCPar.seq,MCMCPar.n);

if strcmp(Update,'Parallel_Direction_Update'), % PARALLEL DIRECTION UPDATE
    % Define which points of Zoff to use to generate jumps
    rr(1,1) = 1; rr(1,2) = rr(1,1) + DEversion(1) - 1; rr(1,3) = rr(1,2) + 1; rr(1,4) = rr(1,3) + DEversion(1) - 1;
    % Do this for each chain
    for qq = 2:MCMCPar.seq,
        % Define rr to be used for population evolution
        rr(qq,1) = rr(qq-1,4) + 1; rr(qq,2) = rr(qq,1) + DEversion(qq,1) - 1; rr(qq,3) = rr(qq,2) + 1; rr(qq,4) = ...
            rr(qq,3) + DEversion(qq,1) - 1;
    end;

    % Each chain evolves using information from other chains to create offspring
    for qq = 1:MCMCPar.seq,

        % ------------ WHICH DIMENSIONS TO UPDATE? USE CROSSOVER ----------
        [i] = find(D(qq,1:MCMCPar.n) > (1-CR(qq,1)));

        % Update at least one dimension
        if isempty(i), i = randperm(MCMCPar.n); i = i(1); end;
        % -----------------------------------------------------------------

        % Generate a random number between 0 and 1
        U = rand;

        % Select the appropriate JumpRate and create a jump
        if (U < MCMCPar.pJumpRate_one),
            % Select the JumpRate (dependent of NrDim and number of pairs)
            NrDim = size(i,2); Gamma = Table_JumpRate(NrDim,DEversion(qq,1));
            % Produce the difference of the pairs used for population evolution
            jump = sum(Zoff(rr(qq,1):rr(qq,2),1:MCMCPar.n) - Zoff(rr(qq,3):rr(qq,4),1:MCMCPar.n),1);
            % Then fill update the dimension
            xjump(qq,i) = Gamma * (1 + noise_x(qq,i)).*jump(1,i);
        else
            % Full space jump
            Gamma = 1; CR(qq,1) = 1;
            % Compute delta from one pair
            jump = Zoff(rr(qq,1),1:MCMCPar.n) - Zoff(rr(qq,4),1:MCMCPar.n);
            % Now jumprate to facilitate jumping from one mode to the other in all dimensions
            xjump(qq,1:MCMCPar.n) = Gamma * jump(1,1:MCMCPar.n);
        end;
    end;
end;

if strcmp(Update,'Snooker_Update'), % SNOOKER UPDATE
    % Determine the number of rows of Zoff
    NZoff = size(Zoff,1);
    % Define rr
    rr = [1:1:NZoff]; rr = reshape(rr,2,size(rr,2)/2)';
    % Define JumpRate -- uniform rand number between 1.2 and 2.2
    Gamma = 1.2 + rand;
    % Loop over the individual chains
    for qq = 1:MCMCPar.seq,
        % Define which points of Zoff z_r1, z_r2
        zR1 = Zoff(rr(qq,1),1:MCMCPar.n); zR2 = Zoff(rr(qq,2),1:MCMCPar.n);
        % Now select z from Zoff; z cannot be zR1 and zR2
        r = [1:NZoff]; r(rr(qq,1)) = 0; r(rr(qq,2)) = 0; r = r(r>0); t = randperm(NZoff-2);
        % Define z
        z(qq,1:MCMCPar.n) = Zoff(r(t(1)),1:MCMCPar.n);
        % Define projection vector x(qq) - z
        F = X(qq,1:MCMCPar.n) - z(qq,1:MCMCPar.n); D = max(F*F',1e-300);
        % Orthogonally project of zR1 and zR2 onto F
        zP = F * (sum((zR1-zR2).*F) / D);
        % And define the jump
        xjump(qq,1:MCMCPar.n) = Gamma * zP;
        % Update CR because we only consider full dimensional updates
        CR(qq,1) = 1;
    end;
end;

% Now propose new x
xnew = X + xjump;

% Define alpha_s
if strcmp(Update,'Snooker_Update'),
    % Determine Euclidean distance
    alpha_s = [(sum((xnew - z).^2,2)./sum((X - z).^2,2)).^((MCMCPar.n-1)/2)];
else
    alpha_s = ones(MCMCPar.seq,1);
end;

% Do boundary handling -- what to do when points fall outside bound
if strcmp(BoundHandling,'Reflect');
    [xnew] = ReflectBounds(xnew,ParRange);
end;

if strcmp(BoundHandling,'Bound');
    [xnew] = SetToBounds(xnew,ParRange);
end;

if strcmp(BoundHandling,'Fold');
    [xnew] = FoldBounds(xnew,ParRange);
end;

if strcmp(BoundHandling,'None');
    % Do nothing
end;