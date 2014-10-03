function MainLabColumn

    close all

    addpath('../../../Common/');
    addpath('../Data/');
    
%     RESULTS_FILE_NAME_TEMPLATE = '../Data/ShirishExperiment_%s.mat';
%     DATA_FILE_NAME = '../../../Investigations/Transforms/mat/ShirishExperiment.mat';
%     mu = 4.54;
%     sigma = 0.775;
%     delay = 400;
%     cAdj = -3.05;
%     tAdj = 1;
%     pvTot = 1.3;
%     kExch = 1.3e-3;
%     beta = Inf;
%     cAdj = 3.5;
    RESULTS_FILE_NAME_TEMPLATE = '../Data/ShirishSimulation_%s.mat';
    DATA_FILE_NAME = '../../../Investigations/Transforms/mat/ShirishSimulation.mat';
    mu = log(132);
    sigma = 0.62;
    delay = 0;
    tAdj = 2;
    pvTot = 6e+0;            % Total pore volume
    kExch = 5e-3;
    beta = 2;               % Immobile-mobile volume ratio
    cAdj = 0;
    fluxAdj = 1e+2;
    
    Const = DefineConstants();

    % Dimensions
    zTop = 0;
    zBottom = -1;
    dz = 0.05;
    ModelDim = InitializeNodes('z', zBottom, zTop, dz);
    
    % Load all data:
    
    % Load: rainData, TimeParams, StartDate
    RawData = load(DATA_FILE_NAME);
    RawData.fluxIn = RawData.fluxIn * fluxAdj;
    RawData.fluxOut = RawData.fluxOut * fluxAdj;
    RawData.t = RawData.t * tAdj;
    
%     % Try to estimate parameters mu and sigma:
%     paramsRange = struct();
%     paramsRange.rangeMu = [4, 5];
%     paramsRange.deltaMu = 0.01;
%     paramsRange.rangeSigma = [0.76, 0.86];
%     paramsRange.deltaSigma = 0.005;
%     [muEst, sigmaEst, err] = LogNormalParamsFromFlux(RawData.t', ...
%         RawData.fluxIn', RawData.fluxOut', paramsRange);
%     % mu = 4.54; sigma = 0.775; err = 0.053798443378235;
%     dt = RawData.t(2) - RawData.t(1);
%     fluxOutNum = NumericalConvolution(RawData.t', RawData.fluxIn', ...
%         @(t, mu_, sigma_) lognpdfX(t, mu_, sigma_, dt), muEst, sigmaEst);
%     plot(RawData.t', [RawData.fluxOut', fluxOutNum]);
    
    % Precipitation in meters per time interval dt
    rainData = 1e-2 * RawData.fluxIn;
%     rainData(1:end) = 0;
%     rainData([10:80, 80:100]) = 1e-2;
%     rainData([6:10,30:35]) = 1e-2;
    % Structure containing time parameters
    TimeParams = TimeParamsCl(RawData.t);
    
    % Concentration of solutes in rainwater
    rainConcentrationData = 0 * ones(size(rainData));
    
    RainInfo = struct();
    RainInfo.intensity = rainData;
    RainInfo.concentration = rainConcentrationData;
    
    %% Model parameters
    ModelParams = struct();
    ModelParams.DO_BIOCHEMISTRY = false;
    ModelParams.DO_RECIRCULATION = false;
    % Parameters of log-normally distributed flow response
    LogNorm = LogNormalCl(mu, sigma, delay, Const);
    ModelParams.LogNorm = LogNorm;
    % Pore volume
    ModelParams.totalPv = pvTot;
    % Immobile-mobile volume ratio
    ModelParams.beta = beta;
    % Source/sink rate
    ModelParams.lambda = 0;
    % Exchange rate between mobile-immobile phases
    ModelParams.kExch = kExch;
    % Initial concentration of intert specie(s)
    ModelParams.mInertIni = 2.3 * [1.23 / 3.076923076923125];

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
        ModelOutput = Compute(TimeParams, RainInfo, ModelDim, ModelParams);
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
    iSpecies = 24;
    ShowPlots(ModelOutput, ModelParams, TimeParams, iSpecies);
    
    % Add experimental data to plots
    % out flux
    PlotAdd(figure(2), RawData.t, RawData.fluxOut, 'Experimental flux');
    % out concentration
    PlotAdd(figure(4), RawData.t, ModelOutput.cRemaining(1, 1, iSpecies) / RawData.c(1) * ...
        (RawData.c + cAdj), 'Experimental result');
    
    fH = figure(5);
    colorOrder = [ ...
        0.6, 0.6, 0.6; ...
        0.0, 0.7, 0.0; ...
        0.3, 0.0, 0.0];
    set(gca, 'ColorOrder', colorOrder, 'NextPlot', 'replacechildren');
    plotyy(RawData.t, [RawData.fluxIn; RawData.fluxOut], RawData.t, RawData.c);
    axH = get(fH, 'children');
    lH = get(axH, 'children');
    range = [min(RawData.fluxIn), max(RawData.fluxIn)];
    rangeWidth = range(2) - range(1);
    range = [range(1) - 0.1 * rangeWidth, range(2) + 0.1 * rangeWidth];
    ylim(range);
    ylabel('Water flux');
    set(get(axH(1), 'ylabel'), 'string', 'Out concentration');
    xlabel('time');
    legend('In flux', 'Out flux', 'Out concentration');
    set(axH(1), 'xlim', [0, 9000]);
    set(axH(2), 'xlim', [0, 9000]);
    
    
    return
    
    close all
    iSpecies = 22; % [1:12, 17:20]; % [2:5, 8:9, 22] % [2, 4]
    PLOT_SEPARATELY = true;
    AnalyzeOutConcentrations(ModelOutput, TimeParams, ModelParams, ...
        ParameterOfInterest, iSpecies, Const, PLOT_SEPARATELY, TimeParams.maxDays);
    
    return
    
    function pdf = lognpdfX(t, mu, sigma, dt)
        pdf = logncdf([t; t(end)+dt], mu, sigma);
        pdf = diff(pdf);
    end
end