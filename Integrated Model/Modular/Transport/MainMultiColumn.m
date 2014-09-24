function MainMultiColumn

    close all

    addpath('../../../Common/');
    addpath('../Data/');
    
    RESULTS_FILE_NAME_TEMPLATE = '../Data/baseline_multi_%s.mat';
    
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
    rainData(1:5) = 1e-2;
    rainData(6:end) = 0;
    
    % Structure containing time parameters
    TimeParams = PrecipitationData.TimeParams;

    % %% For recirculation / irrigation
    % rainData = zeros(size(rainData));
    % totalRecirculation = 3;
    % recirculationIntervalsPerDay = 4;
    % nRecirculationIntervals = TimeParams.maxDays * recirculationIntervalsPerDay;
    % recirculationRate = totalRecirculation / nRecirculationIntervals;
    % iRecirculation = mod(1:TimeParams.numIntervals, 24) < recirculationIntervalsPerDay;
    % rainData(iRecirculation) = recirculationRate;
    % %%

    % Copy rain inputs for one more year
    rainData = 1e+2 * cat(2, rainData, rainData);
    % Update time parameters for an extended interval
    TimeParams.numIntervals = numel(rainData);
    TimeParams.maxDays = TimeParams.numIntervals / TimeParams.intervalsPerDay;
    TimeParams.t = (1:TimeParams.numIntervals) * TimeParams.dt;
    TimeParams.daysElapsed = TimeParams.t / (TimeParams.dt * TimeParams.intervalsPerDay);

    % Concentration of solutes in rainwater
    rainConcentrationData = 0 * ones(size(rainData));
    
    %% Model parameters
    nCells = 2;
    ModelParams = cell(nCells, 1);
    ModelOutput = cell(nCells, 1);
    
    muAll = [0, 1];
    for iCell = 1:nCells
        ModelParams{iCell} = struct();
        ModelParams{iCell}.DO_BIOCHEMISTRY = true;
        ModelParams{iCell}.DO_RECIRCULATION = false;
        mu = muAll(iCell);
        sigma = 1;
        delay = 0;
        LogNorm = LogNormalCl(mu, sigma, delay, Const);
        ModelParams{iCell}.LogNorm = LogNorm;
        % Pore volume
        ModelParams{iCell}.totalPv = 1;
        % Immobile-mobile volume ratio
        ModelParams{iCell}.beta = Inf;
        % Source/sink rate
        ModelParams{iCell}.lambda = 0 * 1e-4;
        % Exchange rate between mobile-immobile phases
        ModelParams{iCell}.kExch = 1e-3;
        % Initial concentration of intert specie(s)
        ModelParams{iCell}.mInertIni = 1;
    end

    TimeParams.maxDays = 40;

    ParameterOfInterest = struct();
    % 'baseline'; 'no_rain'; 'no_chem'; 'short'; 'constant_rain'; 'recirculation'
    ParameterOfInterest.name = 'baseline';
    
    %% Sensitivity analysis of parameters (Name as in code (this name will be also in file name))
%     ModelParams.mu = -1;
%     ParameterOfInterest.name = 'mu=-1';
%     ParameterOfInterest.name = 'beta';
%     ParameterOfInterest.name = 'kExch';
%     ParameterOfInterest.name = 'kExchPart';
%     ParameterOfInterest.name = 'recirculation_clean_water';

    RESULTS_FILE_NAME = sprintf(RESULTS_FILE_NAME_TEMPLATE, ...
        GenerateCharacteristicSuffix(ModelParams{1}, ParameterOfInterest));

    % %% Some test cases
    % % Case #1:
    % %   initial concentrations = [0; 0]
    % %   constant rain, injection of solute at initial time steps
    % cIni = [0; 0];
    % rainData = 1e-2 * ones(size(rainData));
    % rainConcentrationData(1:5) = 1;
    % Case #2:
    %   initial concentrations = [1; 1]
    %   short clean rain at initial time steps
    % rainData(:) = 1e-3;
    % rainData(:) = 0;
    % rainData(1:5) = 1e-3;
    % rainData(361:365) = 1e-3;
    % %% End test cases

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
        for iCell = 1:nCells
            fprintf('> Cell #%d\n', iCell);
            ModelOutput{iCell} = ComputeTravelTimes(TimeParams, rainData, rainConcentrationData, ...
                ModelDim, ModelParams{iCell});
        end
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
    
    % Summarize results
    ModelOutput = SummarizeResults(ModelOutput, ModelParams);
    
    % Validate
    COMP_VARS = {'cOutTotal', 'mOutTotal', 'mRemRes'};
    CheckResults(ModelOutput, ModelParams, validateAction, RESULTS_FILE_NAME, COMP_VARS);

%     %% TODO: EC calculation
%     %    [Ionic strength] = 0.5 * Sum(C_i * Charge_i^2);
%     %    EC = 35.69 * [Ionic strength] + 5.45
%     %  or
%     %    Ec = Sum(C_i * [Specific conductivity])
%     %  Main ions are: Ca2+, Na+, Cl-, NH4+, HCO3-, H+, OH-, SO_4_2-, VFA
%     % Calculate electrical conductivity
%     % Ion conductivity table (name, index of specie, conductivity value)
%     ionCondTable = {...
%         'Cl-', , 76.35
%     }
    
    %% Plotting
    iSpecies = 22;
    ShowPlots(ModelOutput, ModelParams, TimeParams, iSpecies);
    
    close all
    iSpecies = 22; % [1:12, 17:20]; % [2:5, 8:9, 22] % [2, 4]
    PLOT_SEPARATELY = true;
    AnalyzeOutConcentrations(ModelOutput, TimeParams, ModelParams, ...
        ParameterOfInterest, iSpecies, Const, PLOT_SEPARATELY, TimeParams.maxDays);
end