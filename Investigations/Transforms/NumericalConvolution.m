function output = NumericalConvolution(tSpan, input, SystemFun, varargin)
    % Function for performing numerical convolution
    
    output = zeros(size(input));
    
    for iT = 1:numel(tSpan)
        t = tSpan(iT);
        tAfter = tSpan(iT:end) - t;
        if (input(iT) ~= 0)
            output(iT:end) = output(iT:end) + input(iT) * SystemFun(tAfter, varargin{:});
        end
    end
end