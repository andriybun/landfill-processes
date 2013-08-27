function Const = DefineConstants()
    Const = struct();
    
    % Tolerance parameters
    Const.EPSILON = 1e-10;
    Const.CONCENTRATION_EPSILON = 1e-6;
    Const.NUM_SIGMAS = 6;

    % Calculate for results or analyze stored results
    Const.CALCULATE_RESULTS = 1100;
    Const.LOAD_SAVED_RESULTS = 1101;
    
    % Validation action types
    Const.NO_VALIDATION = 1200;
    Const.SAVE_RESULTS = 1201;
    Const.COMPARE_RESULTS = 1202;
    
    
end