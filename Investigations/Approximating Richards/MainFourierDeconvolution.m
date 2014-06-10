function MainFourierDeconvolution

    close all

    addpath('../../Common');
    
    MAT_FILE_DIR = '../Transforms/mat/';
%     CASE_NAME = 'Richards_Output_Pulse_Very_Slow_Response_zref_-1.0.mat';
%     CASE_NAME = 'Richards_Output_Very_Slow_Response.mat'; % Sin, Square_Wave
    caseNameVector = {'Extreme', 'Middle_Extreme', 'Netherlands', 'Uniform'};
    caseIdx = 3;
    CASE_NAME = sprintf('Richards_Output_rainfall_%s.mat', caseNameVector{caseIdx});
    
    RawData = load([MAT_FILE_DIR CASE_NAME]);
    
    % Get data for processing
    inpSel =  ':'; % ':'; % 2125:5000 6500:8760
    t = reshape(RawData.t(inpSel), [], 1);
    t = t - t(1);
    dt = RawData.t(2) - RawData.t(1);
    qIn = -sum(RawData.qIn(inpSel, :), 2) / 15;
    qOut = -sum(RawData.qOut(inpSel, :), 2) / 15;
    
    qInCum = [0; cumsum(qIn)];
    qOutCum = [0; cumsum(qOut)];
    nEl = numel(qInCum);
    step = 25;
    iSel = 1:step:nEl;
    t = t(iSel(1:end-1));
    dt = (t(2) - t(1));
    qIn = diff(qInCum(iSel));
    qOut = diff(qOutCum(iSel));
    
    % Fourier transforms of inputs
    qInF = fft(qIn);
    qOutF = fft(qOut);
    
    % Deconvolution
    pdfF = qOutF ./ qInF;
    
    % Inverse Fourier transform to get approximation of probability density function
    pdfEst = smooth(ifft(pdfF), 1);


%     pdfEst = lognpdfX(t, muEst, sigmaEst, dt);
    
    resSel = 1:40; % ':';
%     plot(pdfEst(resSel));
    pdfEstAdj = pdfEst / sum(pdfEst(resSel));
    [muEst, sigmaEst] = CalcLognormMuSigma(t(resSel), pdfEstAdj(resSel));
%     muEst = 9.4;
%     sigmaEst = 0.55;
    fprintf('%f\t%f\n', muEst, sigmaEst);
    
    % Calculate numerical convolutions
    qOutConvLognVar = NumericalConvolution(t, qIn, pdfEstAdj);
    
    qOutConvLogn = NumericalConvolution(t, qIn, ...
        @(t, mu, sigma) lognpdfX(t, mu, sigma, dt), muEst, sigmaEst);
    thetaIni = 0.157630;
    vPore = 0.45 / 3;
%     [qOutConvLognVar, dTheta, dMu, dSigma] = NumericalConvolution(t, qIn, ...
%         @(t, mu, sigma) lognpdfX(t, mu, sigma, dt), @LogNormalParams, thetaIni, vPore);
%     plot(t, dTheta);

%%
    fH = figure(1);
%     subplot(1, 2, 1); 
    set(fH, 'Position', [100, 100, 390, 300]);
    set(gca, 'FontSize', 8);
    lineH = plot(t(resSel), ...
        cat(2, pdfEstAdj(resSel) / dt, lognpdf(t(resSel), muEst, sigmaEst)) / step);
    set(lineH(1), 'LineWidth', 2, 'Color', [0, 0, 0]);
    set(lineH(2), 'LineWidth', 2, 'Color', [0.6, 0.6, 0.6], 'LineStyle', '--');
    legH = legend('PDF estimated by deconvolution', ...
        ['Fitted log-normal PDF' char(10) '\mu = ' sprintf('%3.5f', muEst) ...
        ', \sigma = ' sprintf('%3.5f', sigmaEst)], 'Location', 'Northeast');
    set(legH, 'FontSize', 8);
%     ylim([0, 0.03]);
    xlabel('Time, minutes');
    ylabel('Flux, m/min');
    hgsave(fH, sprintf('./fig/%s_PDF.fig', CASE_NAME));
%%

    % Plot
    fH = figure(2);
    set(fH, 'Position', [500, 400, 700, 350]);
    plotyy(t, qIn, t, cat(2, qOutConvLogn, qOutConvLognVar, qOut));
    legend('Influx', 'Convolution (lognormal)', 'Convolution (lognormal, variable)', ...
        'Real data', 'Location', 'Northwest');
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
    set(axH(2), 'xlim', [t(1), t(end)]); % [3, 5.24] * 1e5);
    set(axH(3), 'xlim', [t(1), t(end)]);
    set(lH{2}(1), 'Color', [0., 0., 0.], 'LineWidth', 1);
    set(lH{2}(2), 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--', 'LineWidth', 2);
    set(lH{2}(3), 'Color', [0.7, 0.7, 0.7], 'LineWidth', 2);
    hgsave(fH, sprintf('./fig/%s_(de)convolution.fig', CASE_NAME));
    
    return

    function pdf = lognpdfX(t, mu, sigma, dt)
        pdf = logncdf([t; t(end)+dt], mu, sigma);
        pdf = diff(pdf);
    end

    function [optMu, optSigma, optErr] = CalcLognormMuSigma(x, pdf)
        DELTA_MU = 2e-2;
        DELTA_SIGMA = 5e-3;
        RANGE_MU = [2, 5];
        RANGE_SIGMA = [0.8, 1.2];
        NUM_SIGMAS_CONFIDENCE = 4;
        ERR_EPSILON = 3e-4;
        f = @(x) 1; % (1 - erf((x - 4.e+4) / 1.e+4)) / 2;
        pdf = pdf .* f(x);
        pdf = max(0, pdf); %  / sum(pdf);
        optErr = Inf;
        optMu = nan;
        optSigma = nan;
        for mu_ = RANGE_MU(1):DELTA_MU:RANGE_MU(2);
            for sigma_ = RANGE_SIGMA(1):DELTA_SIGMA:RANGE_SIGMA(2)
                err = SumSquares(pdf, lognpdfX(x, mu_, sigma_, dt));
                if (err < optErr)
                    optErr = err;
                    optMu = mu_;
                    optSigma = sigma_;
                end
            end
        end
%         close all;
%         plotyy(x, pdf, x, f(x));
%         hold on
%         plot(x, lognpdfX(x, optMu, optSigma, dt), 'r');
%         hold off
    end

    function err_ = SumSquares(v1, v2, k)
        if nargin < 3
            k = ones(size(v1));
        end
        err_ = sum((k .* (v1 - v2)) .^ 2);
    end

%     function [mux, sigmax] = CalcLognormMuSigma(x, pdf)
%             meanEst = sum(pdf .* x);
%             varEst = sum(pdf .* x.^2) - meanEst ^ 2;
%             [mux, sigmax] = CalcLognormMuSigmaFromMoments(meanEst, varEst);
%     end
%     
%     function [mux, sigmax] = CalcLognormMuSigmaFromMoments(evn, varn)
%         sigmax = sqrt(log(varn ./ (evn .* evn) + 1));
%         mux = log(evn) - sigmax .* sigmax ./ 2;
%     end

end