function Main
%% TODO: water must flow even if no rain (two domain?)
%%

    close all

    addpath('../../../Common/');
    addpath('../Data/');
    
    RESULTS_FILE_NAME_TEMPLATE = '../Data/baseline_test_travel_times_%s.mat';
    
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

    % % Copy rain inputs for one more year
    % rainData = cat(2, rainData, rainData);
    % % rainData(:) = 0;
    % Update time parameters for an extended interval
    % TimeParams.numIntervals = numel(rainData);
    % TimeParams.maxDays = TimeParams.numIntervals / TimeParams.intervalsPerDay;
    % TimeParams.t = (1:TimeParams.numIntervals) * TimeParams.dt;
    % TimeParams.daysElapsed = TimeParams.t / (TimeParams.dt * TimeParams.intervalsPerDay);

    % Concentration of solutes in rainwater
    rainConcentrationData = 0 * ones(size(rainData));
    
    %% Model parameters
    ModelParams = struct();
    % Parameters of log-normally distributed flow response
    ModelParams.mu = 0;
    ModelParams.sigma = 1;
    % Pore volume
    ModelParams.totalPv = 0.44;
    % Immobile-mobile volume ratio
    ModelParams.beta = 10;
    % Source/sink rate
    ModelParams.lambda = 0 * 1e-4;
    % Exchange rate between mobile-immobile phases
    ModelParams.kExch = 1e-1;
    % Exchange rate between mobile phase and particles flowing
    ModelParams.kExchPart = (log(1.5) / exp(ModelParams.mu - ModelParams.sigma ^ 2));

    TimeParams.maxDays = 40;

    ParameterOfInterest = struct();
    ParameterOfInterest.name = 'baseline';
    
    %% Sensitivity analysis of parameters (Name as in code (this name will be also in file name))
%     ModelParams.mu = -1;
%     ParameterOfInterest.name = 'mu=-1';
%     ParameterOfInterest.name = 'beta';
%     ParameterOfInterest.name = 'kExch';
%     ParameterOfInterest.name = 'kExchPart';
%     ParameterOfInterest.name = 'recirculation_clean_water';

    RESULTS_FILE_NAME = sprintf(RESULTS_FILE_NAME_TEMPLATE, ...
        GenerateCharacteristicSuffix(ModelParams, ParameterOfInterest));

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
    rainData = 0 * rainData;
    rainData(1:5) = 1e-3;
    rainData(361:365) = 1e-3;
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
    ShowPlots(ModelOutput, rainData, ModelParams, TimeParams, iSpecies);
    
    close all
    for iSpecies = [2:5, 8:9, 22]
        AnalyzeOutConcentrations(ModelOutput, TimeParams, ModelParams, ...
            ParameterOfInterest, iSpecies);
    end
end