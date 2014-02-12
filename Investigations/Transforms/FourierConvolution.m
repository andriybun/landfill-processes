function output = FourierConvolution(tSpan, input, SystemFun, varargin)
    % Function for performing fourier convolution
    
    inputF = fft(input);
    systemFunF = fft(SystemFun(tSpan, varargin{:}));
    outputF = inputF .* systemFunF;
    output = ifft(outputF);
    
end