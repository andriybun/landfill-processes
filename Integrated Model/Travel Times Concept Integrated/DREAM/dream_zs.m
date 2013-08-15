function [Sequences,Reduced_Seq,X,Z,output] = dream_zs(MCMCPar,ParRange,Measurement,ModelName,Extra,option,Restart)
% ------------- DREAM with sampling from past and snooker updates: DREAM_ZS --------------------%
%                                                                                               %
% The code presented herein is a Markov Chain Monte Carlo algorithm that runs multiple chains   %
% in parallel for efficient posterior exploration. The algorithm, entitled DREAM_(ZS) is        %
% based on the original DREAM sampling scheme, but uses sampling from an archive of past        %
% states to generate candidate points in each individual chain. Theory and numerical examples   %
% of DREAM_(ZS) have been presented in Vrugt et al. (2011). Details can also be found in        %
% Ter Braak and Vrugt (2008).                                                                   %
%                                                                                               %
% Sampling from past has three main advantages:                                                 %
% (1) Circumvents the requirement of using N = d for posterior exploration. This will speed-up  %
% convergence to a limiting distribution, especially for high-dimensional problems (large d).   %
% (2) Outlier chains do not need explicit consideration. By sampling historical states,         %
% aberrant trajectories an jump directly to the modal region at any time during the             %
% simulation. The N path ways simulated with DREAM_(ZS) therefore maintain detailed balance at  %
% every singe step in the chain.                                                                %
% (3) The transition kernel defining the jumps in each of the chains does not require           %
% information about the current states of the chains. This is of great advantage in a           %
% multi-processor environment where the N candidate points can be generated simultaneously so   %
% that each chain can evolve most efficiently on a different computer. Details of this will be  %
% given in a later publication, which should be ready within the next few months.               %
%                                                                                               %
% DREAM_(ZS) also contains a snooker updater to maximize the diversity of candidate points      %
% and generate jumps beyond parallel direction updates. Finally, DREAM_(ZS) contains subspace   %
% learning in a similar way as DREAM, to maximize the squared jumping distance between two      %
% subsequent points in each chain. This idea has been presented in Vrugt et al. (2008) and      %
% shown to significantly increase the efficiency of posterior exploration. All these options    %
% can be activated from the input file.                                                         %
%                                                                                               %
% DREAM_(ZS) developed by Jasper A. Vrugt and Cajo ter Braak                                    %
%                                                                                               %
% This algorithm has been described in:                                                         %
%                                                                                               %
%   C.J.F. ter Braak, and J.A. Vrugt, Differential Evolution Markov Chain with snooker updater  %
%       and fewer chains, Statistics and Computing, 10.1007/s11222-008-9104-9, 2008             %
%                                                                                               %
%   Vrugt, J.A., E. Laloy, and C.J.F. ter Braak, DiffeRential Evolution Adaptive Metropolis     %
%       with Sampling from the Past and Subspace Updating, SIAM journal on Optimization,        %
%       In Review                                                                               %
%                                                                                               %
% For more information please read:                                                             %
%                                                                                               %
%   Vrugt J.A., H.V. Gupta, W. Bouten and S. Sorooshian, A Shuffled Complex Evolution           %
%       Metropolis algorithm for optimization and uncertainty assessment of hydrologic model    %
%       parameters, Water Resour. Res., 39 (8), 1201, doi:10.1029/2002WR001642, 2003.           %
%                                                                                               %
%   ter Braak, C.J.F., A Markov Chain Monte Carlo version of the genetic algorithm Differential %
%       Evolution: easy Bayesian computing for real parameter spaces, Stat. Comput., 16,        %
%       239 - 249, doi:10.1007/s11222-006-8769-1, 2006.                                         %
%                                                                                               %
%   Vrugt, J.A., C.J.F. ter Braak, M.P. Clark, J.M. Hyman, and B.A. Robinson, Treatment of      %
%       input uncertainty in hydrologic modeling: Doing hydrology backward using Markov         %
%       chain Monte Carlo, Water Resour. Res., 44, W00B09, doi:10.1029/2007WR006720, 2008.      %
%                                                                                               %
%   Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon, B.A. Robinson, and J.M. Hyman,       %
%       Accelerating Markov chain Monte Carlo simulation by self adaptive differential          %
%       evolution with randomized subspace sampling, International Journal of Nonlinear         %
%       Sciences and Numerical Simulation, 10(3), 273-290.                                      %
%                                                                                               %
% Copyright (c) 2008, Los Alamos National Security, LLC                                         %
%                                                                                               %
% All rights reserved.                                                                          %
%                                                                                               %
% Copyright 2008. Los Alamos National Security, LLC. This software was produced under U.S.      %
% Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is     %
% operated by Los Alamos National Security, LLC for the U.S. Department of Energy. The U.S.     %
% Government has rights to use, reproduce, and distribute this software.                        %
%                                                                                               %
% NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES A NY WARRANTY, EXPRESS OR  %
% IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE. If software is modified to    %
% produce derivative works, such modified software should be clearly marked, so as not to       %
% confuse it with the version available from LANL.                                              %
%                                                                                               %
% Additionally, redistribution and use in source and binary forms, with or without              %
% modification, are permitted provided that the following conditions are met:                   %
% � Redistributions of source code must retain the above copyright notice, this list of         %
%   conditions and the following disclaimer.                                                    %
% � Redistributions in binary form must reproduce the above copyright notice, this list of      %
%   conditions and the following disclaimer in the documentation and/or other materials         %
%   provided with the distribution.                                                             %
% � Neither the name of Los Alamos National Security, LLC, Los Alamos National Laboratory, LANL %
%   the U.S. Government, nor the names of its contributors may be used to endorse or promote    s%
%   products derived from this software without specific prior written permission.              %
%                                                                                               %
% THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND   %
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES      %
% OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS %
% ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, %
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF   %
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)        %
% HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT %
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,       %
% EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                                            %
%                                                                                               %
%                                                                                               %
% Copyright (c) 2008, Los Alamos National Security, LLC                                         %
%                                                                                               %
% Written by Jasper A. Vrugt: vrugt@lanl.gov                                                    %
%                                                                                               %
% Version 0.5: January 2009                                                                     %
% Version 1.0: April 2011         Maintenance update, explicit treatment of prior distribution  %
%                                                                                               %
% --------------------------------------------------------------------------------------------- %

if strcmp(Restart,'No'),

    % Set random number generator
    opts.Seed = 'sum(100*clock)  % evaluated if it is a string';
    % Then generate new seed
    if ischar(opts.Seed)
        randn('state', eval(opts.Seed));     % random number generator state
    else
        randn('state', opts.Seed);
    end

    % Step 0: Initialize variables
    [MCMCPar,pCR,lCR,CR,Iter,teller,new_teller,Sequences,Z,Table_JumpRate,Reduced_Seq,iloc,output] = ...
        InitVariables(MCMCPar,Extra);
    
    newPb = progressbar(MCMCPar.ndraw, 'Iterations');
    
    % Step 1: Sample MCMCPar.m0 points in the parameter space and store in Z
    if strcmp(Extra.InitPopulation,'LHS_BASED'),
        % Latin hypercube sampling
        [Zinit] = LHSU(ParRange.minn,ParRange.maxn,MCMCPar.m0 + MCMCPar.seq);
    elseif strcmp(Extra.InitPopulation,'COV_BASED');
        % Covariance based sampling
        [Zinit] = repmat(Extra.muX,MCMCPar.m0 + MCMCPar.seq,1) + randn(MCMCPar.m0 + MCMCPar.seq,MCMCPar.n) * chol(Extra.qcov);
        % Do boundary handling -- all points need to be in bound
    elseif strcmp(Extra.InitPopulation,'PRIOR');
        % Create the initial position of each chain by drawing each parameter individually from the prior
        for qq = 1:MCMCPar.n,
            for zz = 1:MCMCPar.m0 + MCMCPar.seq,
                Zinit(zz,qq) = eval(char(Extra.prior(qq)));
            end;
        end;
    end;

    % Boundary handling -- make sure points remain in bound
    if strcmp(Extra.BoundHandling,'Reflect');
        [Zinit] = ReflectBounds(Zinit,ParRange);
    end;
    if strcmp(Extra.BoundHandling,'Bound');
        [Zinit] = SetToBounds(Zinit,ParRange);
    end;
    if strcmp(Extra.BoundHandling,'Fold');
        [Zinit] = FoldBounds(Zinit,ParRange);
    end;
    if strcmp(Extra.BoundHandling,'None');
        % Do nothing
    end;

    % Define initial MCMCPar.m0 rows of Z to be initial sample -- posterior density is not needed and thus not evaluated!!
    Z(1:MCMCPar.m0,1:MCMCPar.n) = Zinit(1:MCMCPar.m0,1:MCMCPar.n);

    % Define initial population from last MCMCPar.seq samples of Zinit
    X = Zinit(MCMCPar.m0 + 1:MCMCPar.m0+MCMCPar.seq,1:MCMCPar.n); clear Zinit;

    % Calculate posterior density associated with each value of X
    [p,log_p] = CompDensity(X,MCMCPar,Measurement,ModelName,Extra,option);
    
    % Append X with information about posterior density (or transformation thereof)
    X = [X p(:,1) log_p];

    % Initialize the sequences
    if strcmp(Extra.save_in_memory,'Yes');
        [Sequences] = InitSequences(X,Sequences,MCMCPar);
    else
        [Sequences] = InitSequences(X,[],MCMCPar);
    end;

    if strcmp(Extra.reduced_sample_collection,'Yes');
        % Reduced sample collection
        iloc_2 = 0;
    else
        % Do nothing
    end;

    % Save N_CR in memory and initialize sum_p2
    output.CR(1,1:MCMCPar.nCR+1) = [Iter pCR]; delta_tot = zeros(1,MCMCPar.nCR);

    % Compute the R-statistic of Gelman and Rubin
    [output.R_stat(1,1:MCMCPar.n+1)] = [Iter Gelman(Sequences(1:iloc,1:MCMCPar.n,1:MCMCPar.seq),MCMCPar)];

else % If a restart run is being done: just load the output from the previous ongoing trial

    ORI = Extra.ORI;
    load Optimization.mat;
    
    if Iter >= MCMCPar.ndraw
       MCMCPar.ndraw = 2 * MCMCPar.ndraw;
    end
    
    newPb = progressbar(MCMCPar.ndraw, 'Iterations');
    Extra.ORI = ORI;
    
end;

% Move prior population to posterior population ...
while (Iter < MCMCPar.ndraw),
    tic
    % Initialize totaccept;
    totaccept = 0;

    % Loop a number of times before calculating convergence diagnostic, etc.
    for gen_number = 1:MCMCPar.steps,

        % Initialize teller
        new_teller = new_teller + 1;

        % Define the current locations and associated posterior densities
        [xold,p_xold,log_p_xold] = GetLocation(X,MCMCPar);

        if (MCMCPar.m < 2 * MCMCPar.DEpairs * MCMCPar.seq),
            % The number of elements of Z is not sufficient
            error('size of Z not sufficient to generate offspring with selected MCMCPar.m0, MCMCPar.seq, and MCMCPar.DEpairs');
        else
            % Without replacement draw rows from Z for proposal creation
            R = randsample(MCMCPar.m, 2 * MCMCPar.DEpairs * MCMCPar.seq); Zoff = Z(R,1:MCMCPar.n);
        end;

        % First generate a random number between 0 and 1
        U = rand;
        % Determine to do parallel direction or snooker update
        if (U <= MCMCPar.parallelUpdate),
            Update = 'Parallel_Direction_Update';
        else
            Update = 'Snooker_Update';
        end;

        % Generate candidate points (proposal) in each chain using either snooker or parallel direction update
        [xnew,CR(:,gen_number),alpha_s] = offde(X(1:MCMCPar.seq,1:MCMCPar.n),Zoff,CR(:,gen_number),MCMCPar,Update,...
            Table_JumpRate,ParRange,Extra.BoundHandling);

        % Compute the likelihood of each proposal in each chain
        [p_xnew,log_p_xnew] = CompDensity(xnew,MCMCPar,Measurement,ModelName,Extra,option); p_xnew = p_xnew(:,1);

        % Apply the acceptance/rejectance rule
        [newgen,accept] = metrop(xnew,p_xnew,log_p_xnew,xold,p_xold,log_p_xold,alpha_s,Measurement,MCMCPar,Extra,option);

        % Check whether to add to sequence or to only store current point
        if strcmp(Extra.save_in_memory,'Yes');
            % Define idx based on iloc
            iloc = iloc + 1;
        else
            % Do nothing -- remain at first location
        end

        % Update the location of the chains
        Sequences(iloc,1:MCMCPar.n+2,1:MCMCPar.seq) = reshape(newgen',1,MCMCPar.n+2,MCMCPar.seq);

        % Check whether to store a reduced sample
        if strcmp(Extra.reduced_sample_collection,'Yes');
            if (new_teller == Extra.T),
                % Update iloc_2 and new_teller
                iloc_2 = iloc_2 + 1; new_teller = 0;
                % Reduced sample collection
                Reduced_Seq(iloc_2,1:MCMCPar.n+2,1:MCMCPar.seq) = reshape(newgen',1,MCMCPar.n+2,MCMCPar.seq);
            else
                % Do nothing -- not yet at Extra.T
            end;
        end;

        % Update X
        X = newgen; clear newgen;

        % Compute squared jumping distance for each CR value
        if strcmp(Extra.pCR,'Update');
            % Calculate the standard deviation of each dimension of X
            r = repmat(std(X(1:MCMCPar.seq,1:MCMCPar.n)),MCMCPar.seq,1);
            % Compute the Euclidean distance between new X and old X
            delta_normX = sum(((xold(1:end,1:MCMCPar.n) - X(1:end,1:MCMCPar.n))./r).^2,2);
            % Use this information to update sum_p2 to update N_CR
            [delta_tot] = CalcDelta(MCMCPar,delta_tot,delta_normX,CR(1:MCMCPar.seq,gen_number));
        end;

        % Check whether to append X to Z
        if (mod(gen_number,MCMCPar.k) == 0),
            % Append X to Z
            Z(MCMCPar.m + 1 : MCMCPar.m + MCMCPar.seq,1:MCMCPar.n+2) = X(1:MCMCPar.seq,1:MCMCPar.n+2);
            % Update MCMCPar.m
            MCMCPar.m = MCMCPar.m + MCMCPar.seq;
        end;

        % How many candidate points have been accepted -- for Acceptance Rate
        totaccept = totaccept + sum(accept);

        % Update Iteration
        Iter = Iter + MCMCPar.seq;
        newPb.update(Iter);
    end;

    % Reduce MCMCPar.steps to get rounded iteration numbers
    if (teller == 2), MCMCPar.steps = MCMCPar.steps + 1; end;

    % Store Important Diagnostic information -- Acceptance Rate
    output.AR(teller,1:2) = [Iter 100 * totaccept/(MCMCPar.steps * MCMCPar.seq)];

    % Store Important Diagnostic information -- Probability of individual crossover values
    output.CR(teller,1:MCMCPar.nCR+1) = [Iter pCR];

    % Check whether to update individual pCR values
    if (Iter <= 0.1 * MCMCPar.ndraw);
        if strcmp(Extra.pCR,'Update');
            % Update pCR values
            [pCR] = AdaptpCR(MCMCPar,delta_tot,lCR);
        end;
    end;

    % Generate CR values based on current pCR values
    [CR,lCRnew] = GenCR(MCMCPar,pCR); lCR = lCR + lCRnew;

    % Calculate Gelman and Rubin Convergence Diagnostic
    if strcmp(Extra.save_in_memory,'Yes');
        start_loc = max(1,floor(0.5*iloc)); end_loc = iloc;
        % Compute the R-statistic using 50% burn-in from Sequences
        [output.R_stat(teller,1:MCMCPar.n+1)] = [Iter Gelman(Sequences(start_loc:end_loc,1:MCMCPar.n,1:MCMCPar.seq),...
            MCMCPar)];
    elseif strcmp(Extra.save_in_memory,'No') & strcmp(Extra.reduced_sample_collection,'Yes');
        start_loc = max(1,floor(0.5*iloc_2)); end_loc = iloc_2;
        % Compute the R-statistic using 50% burn-in from Reduced_Seq
        [output.R_stat(teller,1:MCMCPar.n+1)] = [Iter Gelman(Reduced_Seq(start_loc:end_loc,1:MCMCPar.n,1:MCMCPar.seq),...
            MCMCPar)];
    else
        % Do nothing -- no R-statistics is being computed
    end;

    % Update the teller
    teller = teller + 1;

    % Save statement
    varsExceptOri = setdiff(who, 'ORI');
    save('Optimization.mat', varsExceptOri{:});

    % Print Iteration to screen
    Iter
    toc
end;

% Postprocess output from DREAM before returning arguments
[i] = find(sum(Z,2)==0);
% Remove extra rows from Sequences
if isempty(i) == 0,
    % Derive i
    i = i(1) - 1; Z = Z(1:i,1:MCMCPar.n);
end;

[i] = find(sum(Sequences(:,1:end,1),2)==0);
% Remove extra rows from Sequences
if isempty(i) == 0,
    % Derive i
    i = i(1) - 1; Sequences = Sequences(1:i,:,:);
end;

[i] = find(sum(Reduced_Seq(:,1:end),2)==0);
% Remove extra rows from Sequences
if isempty(i) == 0,
    % Derive i
    i = i(1) - 1; Reduced_Seq = Reduced_Seq(1:i,:,:);
end;

[i] = find(sum(output.AR(:,1:end),2)==0);
if isempty(i) == 0,
    % Derive i
    i = i(1) - 1; output.AR = output.AR(1:i,:); output.CR = output.CR(1:i,:);
end;

[i] = find(sum(output.CR(:,1:end),2)==0);
if isempty(i) == 0,
    % Derive i
    i = i(1) - 1; output.CR = output.CR(1:i,:);
end;

newPb.delete();

end