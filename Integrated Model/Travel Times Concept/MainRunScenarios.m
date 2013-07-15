function MainRunScenarios

    close all

    addpath('../../Common/');
    addpath('../Data/');
    
    global NO_VALIDATION SAVE_RESULTS COMPARE_RESULTS
    DefineGlobalVariables();

    % Tolerance parameters
    global EPSILON NUM_SIGMAS 
    EPSILON = 1e-10;
    NUM_SIGMAS = 6;
    
    % Dimensions
    zTop = 0;
    zBottom = -1;
    dz = 0.05;
    ModelDim = InitializeNodes('z', zBottom, zTop, dz);
    
    % Load: rainData, TimeParams, StartDate
    PrecipitationData = load('precipitation');
    % Precipitation in meters per time interval dt
    rainData = PrecipitationData.rainData;
    rainConcentrationData = 0 * ones(size(rainData));
    TimeParams = PrecipitationData.TimeParams;
    
    % Initial concentration (immobile, mobile phases)
    cIni = [1; 1];
    
    %% Model parameters
    ModelParams = struct();
    % Parameters of log-normally distributed flow response
    ModelParams.mu = 0;
    ModelParams.sigma = 1;
    % Pore volume
    ModelParams.totalPv = 1;
    % Immobile-mobile volume ratio
    ModelParams.beta = 10;
    % Source/sink rate
    ModelParams.lambda = 0 * 1e-4;
    % Exchange rate between mobile-immobile phases
    ModelParams.kExch = 1e-2;
    % Exchange rate between mobile phase and particles flowing
    ModelParams.kExchPart = -(log(1.5) / exp(ModelParams.mu - ModelParams.sigma ^ 2));

    %% Parameter to investigate for sensitivity
    Parameter = struct();
%     Parameter.Name = 'beta';
%     Parameter.Values = [0.1, 1, 5, 10];
%     Parameter.Name = 'kExchPart';
%     Parameter.Values = [-0.11, -1.1022, -2.2];
%     Parameter.Name = 'kExch';
%     Parameter.Values = [1e-4, 1e-2, 1e-1];
    Parameter.Name = 'lambda';
    Parameter.Values = [0, 1e-3, 1e-2];
    %%
    
%     TimeParams.maxDays = 30;

    tic
    
    for iVal = 1:numel(Parameter.Values)
        ModelParams.(Parameter.Name) = Parameter.Values(iVal);
        ModelOutput = ComputeTravelTimes(TimeParams, rainData, rainConcentrationData, ...
            ModelDim, ModelParams, cIni);
        
        FILENAME_TEMPLATE = '../Data/results_%s_%5.4f.mat';
        COMP_VARS = {'cOutTotal', 'mOutTotal', 'mRemRes'};
        
        action = SAVE_RESULTS;
        
        resultsFileName = sprintf(FILENAME_TEMPLATE, Parameter.Name, Parameter.Values(iVal));
        CheckResults(ModelOutput, action, resultsFileName, COMP_VARS);
    end

    toc
    
end