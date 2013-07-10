function FourierDeconvolution
    addpath('../../Integrated Model/Data/');

    InputData = load('precipitation');
    % Precipitation in meters per time interval dt
    TimeParams = InputData.TimeParams;
    
    OutputData = load('baseline');
    
    % Get data for processing
    inpSel = 1000:5000;
    t = TimeParams.t(inpSel) - TimeParams.t(inpSel(1));
    dt = TimeParams.dt;
    qIn = 4 * InputData.rainData(inpSel);
    qOut = OutputData.qOutTotal(inpSel);
    
    % Fourier transforms of inputs
    qInF = fft(qIn);
    qOutF = fft(qOut);
    
    % Deconvolution
    pdfF = qOutF ./ qInF;
    
    % Inverse Fourier transform to get approximation of probability density function
    pdfEst = ifft(pdfF);
    pdfEstSmooth = smooth(pdfEst, 1)';
    resSel = 1:200;
    
    [muEst, sigmaEst] = CalcLognormMuSigma(t(resSel), pdfEst(resSel));
    
    % Plot results and compare with real value
    % Actual log-normal parameters
    mu = 0;
    sigma = 1;
    % Analytical PDF
    pdfLogn = lognpdf(t, mu, sigma) * dt;
    plot(t, cat(1, pdfEstSmooth, pdfLogn, lognpdf(t, muEst, sigmaEst) * dt));
    legend('Estimated by deconvolution', 'Analytical', 'Analytical (estimated parameters)');
    NUM_SIGMAS = 2;
    tBounds = LogNormalBounds(mu, sigma, NUM_SIGMAS);
    xlim([0, tBounds(2)]);
    
    
    fprintf('Estimated value of mu =    %f (real %f)\n', muEst, mu);
    fprintf('Estimated value of sigma = %f (real %f)\n', sigmaEst, sigma);
    
    return

    function [mux, sigmax] = CalcLognormMuSigma(x, pdf)
            meanEst = sum(pdf .* x);
            varEst = sum(pdf .* x.^2) - meanEst ^ 2;
            [mux, sigmax] = CalcLognormMuSigmaFromMoments(meanEst, varEst);
    end
    
    function [mux, sigmax] = CalcLognormMuSigmaFromMoments(evn, varn)
        sigmax = sqrt(log(varn ./ (evn .* evn) + 1));
        mux = log(evn) - sigmax .* sigmax ./ 2;
    end

end