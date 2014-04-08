function MainFourierDeconvolution

    close all

    addpath('../../Common');
    
    MAT_FILE_DIR = '../Transforms/mat/';
%     CASE_NAME = 'Richards_Output_Pulse_Very_Slow_Response_zref_-1.0.mat';
    CASE_NAME = 'Richards_Output_Very_Slow_Response.mat'; % Sin, Square_Wave
    
    RawData = load([MAT_FILE_DIR CASE_NAME]);
    
    % Get data for processing
    inpSel =  ':';
    t = reshape(RawData.t(inpSel), [], 1);
    dt = RawData.t(2) - RawData.t(1);
    qIn = -sum(RawData.qIn(inpSel, :), 2) / 15;
    qOut = -sum(RawData.qOut(inpSel, :), 2) / 15;
    
    qInCum = [0; cumsum(qIn)];
    qOutCum = [0; cumsum(qOut)];
    nEl = numel(qInCum);
    step = 1;
    iSel = 1:step:nEl;
    t = t(iSel(1:end-1));
    dt = (t(2 * step) - t(step));
    qIn = diff(qInCum(iSel));
    qOut = diff(qOutCum(iSel));
    
%     % Fourier transforms of inputs
%     qInF = fft(qIn);
%     qOutF = fft(qOut);
%     
%     % Deconvolution
%     pdfF = qOutF ./ qInF;
%     
%     % Inverse Fourier transform to get approximation of probability density function
%     pdfEst = smooth(ifft(pdfF), 9);

    muEst = 9.6543;
    sigmaEst = 0.6;
    pdfEst = lognpdfX(t, muEst, sigmaEst, dt);
    
    resSel = ':';
    pdfEstAdj = pdfEst / sum(pdfEst(resSel));
    [muEst, sigmaEst] = CalcLognormMuSigma(t(resSel), pdfEstAdj(resSel));
    muEst = muEst + 0.;
    sigmaEst = sigmaEst + 0.;
    
    % Calculate numerical convolutions
    qOutConvLogn = NumericalConvolution(t, qIn, ...
        @(t, mu, sigma) lognpdfX(t, mu, sigma, dt), muEst, sigmaEst);
    thetaIni = 0.157630;
    vPore = 0.45 / 3;
    [qOutConvLognVar, dTheta, dMu, dSigma] = NumericalConvolution(t, qIn, ...
        @(t, mu, sigma) lognpdfX(t, mu, sigma, dt), @LogNormalParams, thetaIni, vPore);
    % plot(t, dTheta);

    % Plot
    fH = figure(2);
    set(fH, 'Position', [500, 400, 700, 350]);
    plotyy(t, qIn, t, cat(2, qOut, qOutConvLogn, qOutConvLognVar));
    legend('Influx', 'Real data', 'Convolution (lognormal)', ...
        'Convolution (lognormal, variable)', 'Location', 'Northwest');
    xlabel('Time, minutes');
    ylabel('In flux');
    % Get handles of lines and change their styles
    axH = get(fH, 'children');
    lH = get(axH, 'children');
    qInColor = [0.9, 0.9, 0.9];
    set(lH{3}, 'Color', qInColor);
    set(axH(3), {'ycolor'}, {qInColor / 2});
    set(get(axH(2), 'ylabel'), 'string', 'Out flux');
    set(axH(2), 'ylim', [0, max(cat(1, qOut, qOutConvLogn, qOutConvLognVar)) * 1.05]);
    set(axH(2), 'xlim', [t(1), t(end)]);
    set(axH(3), 'xlim', [t(1), t(end)]);
    set(lH{2}(1), 'Color', [0.3, 0.3, 0.3], 'LineStyle', '--');
    set(lH{2}(2), 'Color', [0.6, 0.6, 0.6]);
    set(lH{2}(3), 'Color', [0., 0., 0.], 'LineWidth', 2);
    hgsave(fH, sprintf('./fig/%s_(de)convolution.fig', CASE_NAME));
    
    return

    function pdf = lognpdfX(t, mu, sigma, dt)
        pdf = logncdf([t; t(end)+dt], mu, sigma);
        pdf = diff(pdf);
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