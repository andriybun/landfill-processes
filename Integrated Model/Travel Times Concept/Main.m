function MainTravelTimes
%% TODO: water must flow even if no rain (two domain?)
%% TODO: mass balance if lambda ~= 0
%%

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
    ModelParams.totalPv = 0.44;
    % Immobile-mobile volume ratio
    ModelParams.beta = 10;
    % Source/sink rate
    ModelParams.lambda = 0 * 1e-4;
    % Exchange rate between mobile-immobile phases
    ModelParams.kExch = 1e-2;
    % Exchange rate between mobile phase and particles flowing
    ModelParams.kExchPart = -(log(1.5) / exp(ModelParams.mu - ModelParams.sigma ^ 2));

    TimeParams.maxDays = 30;

    %% Some test cases
%     % Case #1:
%     %   initial concentrations = [0; 0]
%     %   constant rain, injection of solute at initial time steps
%     cIni = [0; 0];
%     rainData = 1e-2 * ones(size(rainData));
%     rainConcentrationData(1:5) = 1;
%     % Case #2:
%     %   initial concentrations = [1; 1]
%     %   short clean rain at initial time steps
%     rainData = 0 * rainData;
%     rainData(1:5) = 1e-3;
    %% End test cases

%     profile on
    tic
    
    ModelOutput = ComputeTravelTimes(TimeParams, rainData, rainConcentrationData, ...
        ModelDim, ModelParams, cIni);

    toc
    
%     profile off
%     profile viewer
    
%     plotyy(t, qOutTotal(1:nT), t, mOutTotal(1:nT))
%     return
    
    % Validate
    BASELINE_FILE_NAME = '../Data/baseline';
    COMP_VARS = {'cOutTotal', 'mOutTotal', 'mRemRes'};
%     action = SAVE_RESULTS;
    action = COMPARE_RESULTS;
%     action = NO_VALIDATION;
    
    CheckResults(ModelOutput, action, BASELINE_FILE_NAME, COMP_VARS);

    %% Plotting
%     tShow = (TimeParams.daysElapsed > 150) & (TimeParams.daysElapsed < 250);
%     ShowPlots(qOutTotal, mOutTotal, emissionPotential, rainData, ModelParams.lambda, TimeParams, tShow);
    ShowPlots(ModelOutput, rainData, ModelParams.lambda, TimeParams);
    
end