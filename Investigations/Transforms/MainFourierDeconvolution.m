function MainFourierDeconvolution

    close all

    addpath('../../Common');
    
    MAT_FILE_DIR = 'mat/';
%     CASE_NAME = 'ShirishExperiment';
    CASE_NAME = 'Richards_Output_Very_Slow_Response.mat'; % Square_Wave, Sin, Pulse, _zref_-0.5
    
    RawData = load([MAT_FILE_DIR CASE_NAME]);
    
    % Get data for processing
    inpSel =  ':';
    t = reshape(RawData.t(inpSel), [], 1);
    dt = RawData.t(2) - RawData.t(1);
    qIn = -sum(RawData.qIn(inpSel, :), 2);
    qOut = -sum(RawData.qOut(inpSel, :), 2);
    
    qInCum = [0; cumsum(qIn)];
    qOutCum = [0; cumsum(qOut)];
    nEl = numel(qInCum);
    step = 5;
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
    resSel = 1:125; % ':';
    pdfEst = smooth(ifft(pdfF), 1);
    pdfEstAdj = pdfEst(resSel) / sum(pdfEst(resSel));
    plot(pdfEstAdj);
%     [muEst, sigmaEst] = CalcLognormMuSigma(t(resSel), pdfEstAdj(resSel));
%     muEst = muEst - 0.2;
%     sigmaEst = sigmaEst + 0.05;
    muEst = 9.4;
    sigmaEst =  0.55;
    fprintf('%f\t%f\n', muEst, sigmaEst);
    
    % Analytical PDF
    fH = figure(1);
%     subplot(1, 2, 1); 
    set(fH, 'Position', [100, 100, 390, 300]);
    set(gca, 'FontSize', 8);
    lineH = plot(t(resSel), ...
        cat(2, pdfEstAdj(resSel) / dt, lognpdf(t(resSel), muEst, sigmaEst) / step));
    set(lineH(1), 'LineWidth', 2, 'Color', [0, 0, 0]);
    set(lineH(2), 'LineWidth', 2, 'Color', [0.6, 0.6, 0.6], 'LineStyle', '--');
    legH = legend('PDF estimated by deconvolution', ...
        ['Fitted log-normal PDF' char(10) '\mu = ' sprintf('%3.5f', muEst) ...
        ', \sigma = ' sprintf('%3.5f', sigmaEst)], 'Location', 'Southeast');
    set(legH, 'FontSize', 8);
%     ylim([0, 0.03]);
    xlabel('Time, minutes');
    ylabel('Flux, m/min');
    hgsave(fH, sprintf('./fig/%s_PDF.fig', CASE_NAME));
    
    fprintf('Estimated value of mu =    %f\n', muEst);
    fprintf('Estimated value of sigma = %f\n', sigmaEst);
    
    % Calculate numerical convolutions
    qOutConvLogn = NumericalConvolution(t, qIn, ...
        @(t, mu, sigma) lognpdfX(t, mu, sigma, dt), muEst, sigmaEst);
    qOutConvEst = NumericalConvolution(t, qIn, pdfEst(resSel));
    
    fH = figure(2);
%     subplot(1, 2, 2); 
    set(fH, 'Position', [500, 400, 700, 350]);
    plotyy(t, qIn, t, cat(2, qOut, qOutConvLogn, qOutConvEst));
    legend('Influx', 'Real data', 'Convolution (lognormal)', 'Convolution (estimated PDF)', ...
        'Location', 'Northwest');
    xlabel('Time, minutes');
    ylabel('In flux');
    % Get handles of lines and change their styles
    axH = get(fH, 'children');
    lH = get(axH, 'children');
    qInColor = [0.9, 0.9, 0.9];
    set(lH{3}, 'Color', qInColor);
    set(axH(3), {'ycolor'}, {qInColor / 2});
    set(get(axH(2), 'ylabel'), 'string', 'Out flux');
    set(axH(2), 'ylim', [0, max(cat(1, qOut, qOutConvLogn, qOutConvEst)) * 1.05]);
    set(axH(2), 'xlim', [t(1), t(end)]);
    set(axH(3), 'xlim', [t(1), t(end)]);
    set(lH{2}(1), 'Color', [0.6, 0.6, 0.6]);
    set(lH{2}(2), 'Color', [0.3, 0.3, 0.3], 'LineStyle', '--');
    set(lH{2}(3), 'Color', [0., 0., 0.], 'LineWidth', 2);
    hgsave(fH, sprintf('./fig/%s_(de)convolution.fig', CASE_NAME));
    
    fH = figure(3);
    plot(t, cat(2, cumsum(qOut), cumsum(qOutConvLogn), cumsum(qOutConvEst)));
    legend('Real data', 'Convolution (lognormal)', 'Convolution (estimated PDF)', ...
        'Location', 'Southeast');
    ylim([0, max(max([cumsum(qOut), cumsum(qOutConvLogn), cumsum(qOutConvEst)]))]);
    hgsave(fH, sprintf('./fig/%s_cumsum.fig', CASE_NAME));
    
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