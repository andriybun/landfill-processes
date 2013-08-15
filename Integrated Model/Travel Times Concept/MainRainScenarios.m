function MainRainScenarios

    close all

    addpath('../../Common/');
    addpath('../Data/');
    
    Const = DefineConstants();
        
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

    %% Irrigation scenarios
    [~, nT] = size(rainData);
    nScenarios = 4;
    scenarioInfo = cell(1, nScenarios);
    sumRain = sum(rainData);
    nRainIntervals = sum(abs(rainData) > 0);
    rainInterval = nT / nRainIntervals;
    maxRainIntensity = max(abs(rainData));
    rainDataAll = zeros(nScenarios, nT);
    rainDataAll(1, :) = rainData;
    scenarioInfo{1} = 'Real_Data';
    rainDataAll(2, 1:rainInterval:nT) = sumRain / nRainIntervals;
    scenarioInfo{2} = 'Uniform_Intervals';
    rainDataAll(3, 1:nRainIntervals) = sumRain / nRainIntervals;
    scenarioInfo{3} = 'All_at_Once';
    rainDataAll(4, :) = sumRain / nT;
    scenarioInfo{4} = 'Uniform_Intervals';
    % Append some dry period to buffer all delays at the end
    daysExtra = 30;
    ntExtra = daysExtra * TimeParams.intervalsPerDay;
    rainDataAll = cat(2, rainDataAll, zeros(nScenarios, ntExtra));
    rainConcentrationData = cat(2, rainConcentrationData, zeros(1, ntExtra));
    TimeParams.t = cat(2, TimeParams.t, TimeParams.t(end) + TimeParams.dt .* (1:ntExtra));
    TimeParams.numIntervals = TimeParams.numIntervals + ntExtra;
    TimeParams.maxDays = TimeParams.maxDays + daysExtra;
    TimeParams.daysElapsed = TimeParams.t;
    %%
    
%     TimeParams.maxDays = 30;

    tic
    
    for iVal = 1:nScenarios
        rainData = rainDataAll(iVal, :);
        ModelOutput = ComputeTravelTimes(TimeParams, rainData, rainConcentrationData, ...
            ModelDim, ModelParams, cIni);
        
        FILENAME_TEMPLATE = '../Data/results_%s_%s.mat';
        COMP_VARS = {'cOutTotal', 'mOutTotal', 'mRemRes'};
        
        action = Const.SAVE_RESULTS;
        
        resultsFileName = sprintf(FILENAME_TEMPLATE, 'Rain_Data', scenarioInfo{iVal});
        CheckResults(ModelOutput, action, resultsFileName, COMP_VARS);
    end

    toc
    
end