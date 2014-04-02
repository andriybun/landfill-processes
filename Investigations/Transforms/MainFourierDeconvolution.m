function MainFourierDeconvolution

    close all

    addpath('../../Common');
    
    MAT_FILE_DIR = 'mat/';
%     CASE_NAME = 'ShirishExperiment';
    CASE_NAME = 'Richards_Output_Pulse_Very_Slow_Response'; % Square_Wave, Sin, Pulse
    
    
    RawData = load([MAT_FILE_DIR CASE_NAME]);
    
    % Get data for processing
    inpSel = ':';
    t = RawData.t(inpSel);
    dt = RawData.t(2) - RawData.t(1);
    qIn = -sum(RawData.qIn(inpSel, :), 2);
    qOut = -sum(RawData.qOut(inpSel, :), 2);
    
    qInCum = [0; cumsum(qIn)];
    qOutCum = [0; cumsum(qOut)];
    nEl = numel(qInCum);
    step = 2;
    iSel = 1:step:nEl;
    t = t(iSel(1:end-1));
    dt = (t(2 * step) - t(step));
    qIn = diff(qInCum(iSel));
    qOut = diff(qOutCum(iSel));
    
    figure();
    plot(t, [qIn, qOut]);
    legend('input', 'output');
    
    % Fourier transforms of inputs
    qInF = fft(qIn);
    qOutF = fft(qOut);
    
    % Deconvolution
    pdfF = qOutF ./ qInF;
    
    % Inverse Fourier transform to get approximation of probability density function
    pdfEst = smooth(ifft(pdfF), 9);
    
    plot(pdfEst);
    
    resSel = ':';
    pdfEstAdj = pdfEst / sum(pdfEst(resSel));
    [muEst, sigmaEst] = CalcLognormMuSigma(t(resSel), pdfEstAdj(resSel));
    muEst = 9.6543;
    sigmaEst = 0.6;
%     muEst = 9.3;
%     sigmaEst = 0.5;
    
    % Analytical PDF
    fH = figure(1);
%     subplot(1, 2, 1); 
%     set(fH, 'Position', [100, 100, 350, 270]);
%     set(gca, 'FontSize', 8);
    plot(t(resSel), cat(2, pdfEstAdj(resSel), lognpdf(t(resSel), muEst, sigmaEst) * dt / step));
    lH = legend('PDF estimated by deconvolution', ...
        ['Fitted log-normal PDF' char(10) '\mu = ' sprintf('%3.5f', muEst) ...
        ', \sigma = ' sprintf('%3.5f', sigmaEst)]);
    set(lH, 'FontSize', 8);
%     ylim([0, 0.03]);
    xlabel('Time, minutes');
    ylabel('Flux, dimensionless');
    hgsave(fH, sprintf('./fig/%s_PDF', CASE_NAME));
    
    fprintf('Estimated value of mu =    %f\n', muEst);
    fprintf('Estimated value of sigma = %f\n', sigmaEst);
    
    % Calculate numerical convolutions
    qOutConvLogn = NumericalConvolution(t, qIn, ...
        @(t, mu, sigma) lognpdfX(t, mu, sigma, dt), muEst, sigmaEst);
    qOutConvEst = NumericalConvolution(t, qIn, pdfEst(resSel));
    
    fH = figure(2);
%     subplot(1, 2, 2); 
    set(fH, 'Position', [500, 400, 700, 350]);
    plotyy(t, cat(2, qOut, qOutConvLogn, qOutConvEst), t, qIn);
    legend('Real data', 'Convolution (lognormal)', 'Convolution (estimated PDF)', 'Influx', ...
        'Location', 'Northeast');
    xlabel('Time, minutes');
    ylabel('Flux, dimensionless');
    hgsave(fH, sprintf('./fig/%s_(de)convolution', CASE_NAME));
    
    fH = figure(3);
    plot(t, cat(2, cumsum(qOut), cumsum(qOutConvLogn), cumsum(qOutConvEst)));
    legend('Real data', 'Convolution (lognormal)', 'Convolution (estimated PDF)', ...
        'Location', 'Southeast');
    ylim([0, max(max([cumsum(qOut), cumsum(qOutConvLogn), cumsum(qOutConvEst)]))]);
    hgsave(fH, sprintf('./fig/%s_cumsum', CASE_NAME));
    
    return

    function pdf = lognpdfX(t, mu, sigma, dt)
        pdf = logncdf([t; t(end)+dt], mu, sigma);
        pdf = diff(pdf);
%         pdf = lognpdf(t, mu, sigma) * dt;
    end
    
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