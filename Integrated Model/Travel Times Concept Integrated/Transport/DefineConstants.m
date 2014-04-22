function Const = DefineConstants()
    Const = struct();
    
    % Tolerance parameters
    Const.EPSILON = 1e-10;                  % Used for general comparison of floats
    Const.VOLUME_EPSILON = 1e-9;            % Used for comparison of particle volumes
    Const.CONCENTRATION_EPSILON = 1e-6;     % Used for comparison of solute concentrations
    Const.NUM_SIGMAS = 6;

    % Calculate for results or analyze stored results
    Const.CALCULATE_RESULTS = 1100;
    Const.LOAD_SAVED_RESULTS = 1101;
    
    % Validation action types
    Const.NO_VALIDATION = 1200;
    Const.SAVE_RESULTS = 1201;
    Const.COMPARE_RESULTS = 1202;
    
    
end