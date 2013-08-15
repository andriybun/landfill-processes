clc
clear all
close all

pname = '../../Model';
addpath(genpath(pname));

%---------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------
% SET PARAMETERS DDREAM

MCMCPar.ndraw = 2000;           % Maximum number of function evaluations
MCMCPar.parallelUpdate = 0.9;   % Fraction of parallel direction updates
MCMCPar.pJumpRate_one = 0.80;   % Probability of selecting a jumprate of 1 --> jump between modes
MCMCPar.seq = 3;                % Number of Markov Chains / sequences
MCMCPar.DEpairs = 1;            % Number of chain pairs to generate candidate points
MCMCPar.Gamma = 0;              % Kurtosis parameter Bayesian Inference Scheme
MCMCPar.nCR = 3;                % Number of crossover values used
MCMCPar.k = 10;                 % Thinning parameter for appending X to Z
MCMCPar.eps = 5e-2;             % Perturbation for ergodicity
MCMCPar.steps = 10;             % Number of steps before calculating convergence diagnostics
MCMCPar.n = 13;                 % Dimension of the problem (number of parameters to be estimated)
MCMCPar.m0 = 100 * MCMCPar.n;    % Initial size of Z 

ExtraPar.pCR = 'Update';                   % Adaptive tuning of crossover values
ExtraPar.reduced_sample_collection = 'No'; % Thinned sample collection?
ExtraPar.T = 1000;                         % Every Tth sample is collected
ExtraPar.InitPopulation = 'LHS_BASED';     % What type of initial sampling
ExtraPar.BoundHandling = 'Reflect';        % Define the boundary handling 
ExtraPar.save_in_memory = 'Yes';           % Save in memory or not
ExtraPar.FloorParRange = zeros(1,MCMCPar.n); %these are used in the Run_ADE function to scale the dream parameters between 0 and 10 
ExtraPar.FactorParRange = ones(1,MCMCPar.n); %back to the actual parameter values: RealPar = FloorParRange + FactorParRange*DreamPar
ExtraPar.Option = 3; % Define likelihood function (option 6 implemented by S.Korteland, august 2011)
%---------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------
% SET PARAMETERS MSWS1

% Model name
ModelName = 'integrate_bioreactor';
ExtraPar.modus = 1;

% Initial concentration and parameters 
[ExtraPar.IC, ExtraPar.Comp, ExtraPar.Pm, ExtraPar.S, ExtraPar.Rp, ExtraPar.H] = initialize_ODE(0,0);
ExtraPar.ORI = initialize_ORI(ExtraPar.Comp,0);

% Measurement Data
Meas = load_data(ExtraPar.Pm);
Measurement.MeasData = [Meas.Bio_m;Meas.pH_m;Meas.VFA_m;Meas.NH4_m;Meas.CH4p_m;Meas.CO2p_m];
ExtraPar.Dnames = {'Biogas' 'pH' 'VFAx' 'NH4+NH3' 'pCH4' 'pCO2'};
ExtraPar.tm = {Meas.Bio_tm Meas.pH_tm Meas.VFA_tm Meas.NH4_tm Meas.CH4p_tm Meas.CO2p_tm};
Measurement.N = size(Measurement.MeasData(:),1); % number of measurements

% Parameters to be optimized:
% % Case 1
% ExtraPar.Pnames = {'Cxm0'    'Cxs0'  'Ks'      'k'      'qsmaxm' 'qsmaxs' 'pHi(hyd)' 'pHi(meth)'...
%                     'VFAi'    'NH3i'  'H2Sim'   'H2Sis'  'rNH3'    'kla'};
% Parameters    =   [ 0.0022   0.0022   0.00018   0.108    1.79     22.5     5          6.3 ...
%                      0.0234   1.21e-3  0.0047    0.0025   0.0075    2];
%                  
% ParRange.minn =   [0.01*Parameters(1:2) 0.1*Parameters(3)   0.01*Parameters(4:6) ...
%                    0.15*Parameters(7:8) 0.01*Parameters(9:12) 0.01*Parameters(13) 0.01*Parameters(14)];
% 
% ParRange.maxn =   [100*Parameters(1:2)  10000*Parameters(3) 100*Parameters(4:6)  ...
%                    3*Parameters(7:8)    100*Parameters(9:12) 10*Parameters(13) 100*Parameters(14)];
               
% % Case 2
% ExtraPar.Pnames = {'Cxm0'    'Cxs0'  'rNH3'    'kla'};
% Parameters    =   [ 0.0022   0.0022   0.0075    2];
%                  
% ParRange.minn =   [0.01*Parameters(1:2) 0.01*Parameters(3) 0.01*Parameters(4)];
% ParRange.maxn =   [100*Parameters(1:2)  10*Parameters(3)   100*Parameters(4)];

% % Case 3
% ExtraPar.Pnames = {'Cxm0'    'Cxs0'  'Ks'      'k'      'qsmaxm' 'qsmaxs' 'pHi(hyd)' 'pHi(meth)'...
%                     'VFAi'    'NH3i'  'H2Sim'   'H2Sis'  'rNH3'    'kla'};
% Parameters    =   [ 0.0022   0.0022   0.00018   0.108    1.79     22.5     5          6.3 ...
%                      0.0234   1.21e-3  0.0047    0.0025   0.0075    2];
%                  
% ParRange.minn =   [ 0.0001   0.0001   0.0000018 0.0001   0.01     0.01     2          3 ...
%                      0.0002   0.00001  0.00001   0.00001  0.0001   0.01];
% 
% ParRange.maxn =   [ 0.1      0.1      0.018     10       10       50       8          8 ...
%                      0.5      0.5      0.5       0.5      0.01     100];

% Case 4
ExtraPar.Pnames = {'C(hyd)' 'C(Ace)' 'C(NH3)' 'C(CH4)' 'C(SO4)' 'C(H2S)' 'Cxa'  'Cxm'    'Cxs'   'pH2CO3' 'C(Ca)' 'C(Na)' 'C(Cl)'};

Parameters    =   [ 5.075    0        0.065    0         0.07    0.00001   0     0.007    0.0022  -1.5571  1.07    0.1951  0.0951];
                 
ParRange.minn =   [ 1e-6     1e-6     1e-6     0         1e-6    1e-6      0     0        0       -3       1e-6    1e-6    1e-6  ];
 
ParRange.maxn =   [ 10       1        0.1      0.001     0.5     0.01      0.1   0.1      0.1     -1.54    5       0.5     0.1];
%---------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------
% START DDREAM

Restart='Yes';
[Sequences,Reduced_Seq,X,Z,output] = ...
    dream_zs(MCMCPar,ParRange,Measurement,ModelName,ExtraPar,ExtraPar.Option,Restart);
%---------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------
% PLOT RESULTS

plotDREAM(MCMCPar,Z,ParRange,Sequences,output,ExtraPar,Measurement)
