function MainLabColumn
%% TODO: water must flow even if no rain (two domain?)
%% TODO: beta? not taken into account?
%%

    close all

    addpath('../../../Common/');
    addpath('../Data/');
    
    RESULTS_FILE_NAME_TEMPLATE = '../Data/ShirishExperiment_%s.mat';
    DATA_FILE_NAME = '../../../Investigations/Transforms/mat/ShirishExperiment.mat';
    
    Const = DefineConstants();

    % Dimensions
    zTop = 0;
    zBottom = -1;
    dz = 0.05;
    ModelDim = InitializeNodes('z', zBottom, zTop, dz);
    
    % Load all data:
    
    % Load: rainData, TimeParams, StartDate
    RawData = load(DATA_FILE_NAME);
    % Precipitation in meters per time interval dt
    rainData = 1e+2 * RawData.fluxIn;
    rainData(1:5) = 1e-2;
    rainData(6:end) = 0;
    % Structure containing time parameters
    TimeParams.t = RawData.t;
    TimeParams.dt = diff(TimeParams.t(1:2));
    TimeParams.numIntervals = numel(TimeParams.t);
    TimeParams.maxDays = (TimeParams.t(end) - TimeParams.t(1)) / (60 * 60 * 24);
    % TimeParams =
    %
    %             maxDays: 365
    %     intervalsPerDay: 24
    %                  dt: 0.041666666666667
    %        numIntervals: 8760
    %         daysElapsed: [1x8760 double]
    %                   t: [1x8760 double]
    
    % Concentration of solutes in rainwater
    rainConcentrationData = 0 * ones(size(rainData));
    
    %% Model parameters
    ModelParams = struct();
    ModelParams.DO_BIOCHEMISTRY = false;
    ModelParams.DO_RECIRCULATION = false;
    % Parameters of log-normally distributed flow response
    mu = 6;
    sigma = 4.5e-1;
    delay = 400;
    LogNorm = LogNormalCl(mu, sigma, delay, Const);
    ModelParams.LogNorm = LogNorm;
    % Pore volume
    ModelParams.totalPv = 3;
    % Immobile-mobile volume ratio
    ModelParams.beta = Inf;
    % Source/sink rate
    ModelParams.lambda = 0 * 1e-4;
    % Exchange rate between mobile-immobile phases
    ModelParams.kExch = 1e-2;

%     TimeParams.maxDays = 40;

    ParameterOfInterest = struct();
    ParameterOfInterest.name = 'baseline';
    
    RESULTS_FILE_NAME = sprintf(RESULTS_FILE_NAME_TEMPLATE, ...
        GenerateCharacteristicSuffix(ModelParams, ParameterOfInterest));

    resultSource = Const.CALCULATE_RESULTS;
%     resultSource = Const.LOAD_SAVED_RESULTS;

    validateAction = Const.SAVE_RESULTS;
%     validateAction = Const.COMPARE_RESULTS;
%     validateAction = Const.NO_VALIDATION;

    profilerOn = false;
    
    if (resultSource == Const.CALCULATE_RESULTS)
        if profilerOn
            profile on
        end
        
        tic
        % Main computations are done here
        ModelOutput = ComputeTravelTimes(TimeParams, rainData, rainConcentrationData, ...
            ModelDim, ModelParams);
        toc
        
        if profilerOn
            profile off
            profile viewer
        end
    elseif (resultSource == Const.LOAD_SAVED_RESULTS)
        ModelOutput = load(RESULTS_FILE_NAME);
        ModelParams = ModelOutput.ModelParams;
        validateAction = Const.NO_VALIDATION;
    end

    % Validate
    COMP_VARS = {'cOutTotal', 'mOutTotal', 'mRemRes'};
    CheckResults(ModelOutput, ModelParams, validateAction, RESULTS_FILE_NAME, COMP_VARS);
    
    %% Plotting
    iSpecies = 25;
    ShowPlots(ModelOutput, ModelParams, TimeParams, iSpecies);
    figure(4);
    hold on;
    plot(RawData.t, RawData.c-3.05, 'r');
    
    return
    
    close all
    iSpecies = 22; % [1:12, 17:20]; % [2:5, 8:9, 22] % [2, 4]
    PLOT_SEPARATELY = true;
    AnalyzeOutConcentrations(ModelOutput, TimeParams, ModelParams, ...
        ParameterOfInterest, iSpecies, Const, PLOT_SEPARATELY, TimeParams.maxDays);
    
    return
    
    
end