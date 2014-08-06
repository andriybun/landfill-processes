% [LASTMSG, LASTID] = lastwarn;  warning('off',LASTID);
clear all; parallel = 1; pname = '../MSWS1'; addpath(genpath(pname)); addpath('../DREAM');

% Problem specific parameter settings
MCMCPar.n = 14;                          % Dimension of the problem (number of parameters to be estimated)
MCMCPar.ndraw = 30000;                  % Maximum number of function evaluations
MCMCPar.parallelUpdate = 0.9;           % Fraction of parallel direction updates
MCMCPar.pJumpRate_one = 0.20;           % Probability of selecting a jumprate of 1 --> jump between modes

% Recommended parameter settings
MCMCPar.seq = 3;                        % Number of Markov Chains / sequences
MCMCPar.DEpairs = 1;                    % Number of chain pairs to generate candidate points
MCMCPar.Gamma = 0;                      % Kurtosis parameter Bayesian Inference Scheme
MCMCPar.nCR = 3;                        % Number of crossover values used
MCMCPar.m0 = 10 * MCMCPar.n;            % Initial size of Z
MCMCPar.k = 10;                         % Thinning parameter for appending X to Z
MCMCPar.eps = 5e-2;                     % Perturbation for ergodicity
MCMCPar.steps = 10;                     % Number of steps before calculating convergence diagnostics

% --------------------------------------------------------------------------------------------
Extra.pCR = 'Update';                   % Adaptive tuning of crossover values
% --------------------------------------------------------------------------------------------

% --------------------------------------- Added for reduced sample storage -------------------
Extra.reduced_sample_collection = 'No'; % Thinned sample collection?
Extra.T = 1000;                         % Every Tth sample is collected
% --------------------------------------------------------------------------------------------

% What type of initial sampling
Extra.InitPopulation = 'LHS_BASED';
% Give the parameter ranges (minimum and maximum values)
% ParRange.minn = [ 1e-6 0.01  0.02 0.0004 0.00003 0.000178 0.00027 0.00027 1  0.01 0.1 0.001 0.01 0.01];
% ParRange.maxn = [ 0.26 3.13  5.32 0.08   0.042   0.00026  0.1     0.1     30 3    5   0.1   3    3   ];
ParRange.minni = [ 0.09 0.04  0.024 0.001 0.004  0.007 0.27 0.27  1   0.01 0.1 0.001 0.01 0.01];
ParRange.maxni = [ 0.26 0.47  2.4   0.01  0.104  0.7   19   19    30  3    5   0.1   3    3   ];
ParRange.minn = [1e-2.*ParRange.minni];% ParRange.minni(9:end)];
ParRange.maxn = [1e2.*ParRange.maxni];% ParRange.maxni(9:end)];
% Define the boundary handling
Extra.BoundHandling = 'Reflect';
% Save in memory or not
Extra.save_in_memory = 'Yes';

% Load data & initial conditions & model parameters
[Extra.IM, Extra.Pm, Extra.ORI, Extra.Meas, Extra.Comp] = initialize('Pmatrix.csv'); Extra.sw = 1;

% Define the measured data
Measurement.MeasData = []; Measurement.Sigma = []; Measurement.N = size(Measurement.MeasData,1);
% Define modelName
ModelName = 'integrate_bioreactor';
% Define likelihood function -- Sum of Squared Error
option = 4;

% No restart -- just running for the first time. Restart can be used if run
% is termined while running; Use Restart = 'Yes'; to reinitialize the code
% from the saved files.
Restart = 'Yes';
% Name for results file
Extra.fn_DREAM = 'test5_rangepaperx100_weighted_2'; 
if strcmp(Restart,'No'); try load(Extra.fn_DREAM); disp('filename already exists!'); break; catch; end; end;

if parallel == 0, 
    % Run the distributed DREAM algorithm with sampling from past
    [Sequences,Reduced_Seq,X,Z,output] = dream_zs2(MCMCPar,ParRange,Measurement,ModelName,Extra,option,Restart);
    
elseif parallel == 1,
    % Checking if there are any workers still running
    matlabpoolIsOpen = (matlabpool('size') > 0);
    % If yes, close them
    if matlabpoolIsOpen; matlabpool close; end
    % Start matlabpool
    Extra.numWorkers = 3; matlabpool(Extra.numWorkers);
    
    % Run the distributed DREAM algorithm with sampling from past
    [Sequences,Reduced_Seq,X,Z,output] = dream_zs_par2(MCMCPar,ParRange,Measurement,ModelName,Extra,option,Restart);
    
    %  Stop workers
    matlabpool close
end