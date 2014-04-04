function output = NumericalConvolutionVar(tSpan, input, SystemFun, varargin)
    % Function for performing numerical convolution
    
    % Find the type of probability distribution
    switch class(SystemFun)
        case 'double'
            % Vector
            caseId = 1;
            nEl = numel(SystemFun);
            % Normalize it to produce sum = 1
            SystemFun = SystemFun / sum(SystemFun);
        case 'function_handle'
            % Function handle with parameters stored in varargin
            % Check the type of 4th argument
            switch class(varargin{1})
                case 'double'
                    % Numeric value - these parameters are passed to SystemFun
                    caseId = 2;
                case 'function_handle'
                    % Function - this function is used to compute parameters to be passed to
                    % SystemFun
                    caseId = 3;
                otherwise
                    error('Wrong class of parameters to travel time distribution function.');
            end
        otherwise
            error('Wrong class of travel time distribution.');
    end
    
    output = zeros(size(input));
    
    % Initialization of change in volume
    dV = 0;
    nT = numel(tSpan);
    for iT = 1:nT
        t = tSpan(iT);
        tAfter = tSpan(iT:end) - t;
        if (input(iT) ~= 0)
            dV = dV + input(iT);
            switch caseId
                case 1
                    iUpd = iT:min(nT, iT+nEl-1);
                    output(iUpd) = output(iUpd) + input(iT) * SystemFun(1:min(numel(tAfter), nEl));
                case 2
                    varargpass = varargin;
                    output(iT:end) = output(iT:end) + input(iT) * SystemFun(tAfter, varargpass{:});
                case 3
                    theta = varargin{2};
                    vPore = varargin{3};
                    thetaUpd = theta + dV / vPore;
                    params = varargin{1}(thetaUpd);
                    output(iT:end) = output(iT:end) + input(iT) * SystemFun(tAfter, params{:});
                otherwise
                    error('Error');
            end
        end
        dV = dV - output(iT);
    end
end