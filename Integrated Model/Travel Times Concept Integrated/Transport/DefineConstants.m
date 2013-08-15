function Const = DefineConstants()
    Const = struct();
    
    % Tolerance parameters
    Const.EPSILON = 1e-10;
    Const.NUM_SIGMAS = 6;

    % Validation action types
    Const.NO_VALIDATION = 1200;
    Const.SAVE_RESULTS = 1201;
    Const.COMPARE_RESULTS = 1202;
    
    
end