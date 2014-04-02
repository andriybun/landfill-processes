function output = NumericalConvolution(tSpan, input, SystemFun, varargin)
    % Function for performing numerical convolution
    
    % Find the type of probability distribution
    switch class(SystemFun)
        case 'double'
            % Vector
            isFunction = false;
            nEl = numel(SystemFun);
            % Normalize it to produce sum = 1
            SystemFun = SystemFun / sum(SystemFun);
        case 'function_handle'
            % Function handle with parameters stored in varargin
            isFunction = true;
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
            if isFunction
                varargpass = varargin;
                output(iT:end) = output(iT:end) + input(iT) * SystemFun(tAfter, varargpass{:});
            else
                iUpd = iT:min(nT, iT+nEl-1);
                output(iUpd) = output(iUpd) + input(iT) * SystemFun(1:min(numel(tAfter), nEl));
            end
        end
        dV = dV - output(iT);
    end
end