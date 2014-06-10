function MainFourierDeconvolutionMovingWindow
    % Describe function

    %% PARAMETERS
    % Desired tolerance (adjust empirically)
    ERR_EPSILON = 1.5e-4;
    
    % Time stepping parameters
    COARSE_SCALE = 5;               % Coarsen time stepping by a scale of ...
    STEP = 50;                      % Time step for considering windows
    INITIAL_WINDOW_WIDTH = 4000;    % Number of time steps/intervals per window
    
    % Ranges of log-normal parameters
    RANGE_MU = [2, 5];
    RANGE_SIGMA = [0.5, 1.2];
    DELTA_MU = 2e-2;
    DELTA_SIGMA = 5e-3;
    
    % Do calculations or load saved results
    DO_CALCULATIONS = 1;
    LOAD_RESULTS = 2;
    ACTION = LOAD_RESULTS;

    %% END PARAMETERS

    close all

    addpath('../../Common');
    
    MAT_FILE_DIR = '../Transforms/mat/';
    % CASE_NAME = 'Richards_Output_Pulse_Very_Slow_Response_zref_-1.0.mat';
    caseNameVector = {'Extreme', 'Middle_Extreme', 'Netherlands', 'Uniform'};
    caseIdx = 3;
    CASE_NAME = sprintf('Richards_Output_rainfall_%s.mat', caseNameVector{caseIdx});
    % Name of savefile
    SAVE_FILE = sprintf('DeconvoluteMuSigma_%s', caseNameVector{caseIdx});
    
    RawData = load([MAT_FILE_DIR CASE_NAME]);
    
    % Get data for processing
    inpSel =  ':'; % ':'; % 2125:5000 6500:8760
    t = reshape(RawData.t(inpSel), [], 1);
    t = t - t(1);
    dt = RawData.t(2) - RawData.t(1);
    qIn = -sum(RawData.qIn(inpSel, :), 2) / 15;
    qOut = -sum(RawData.qOut(inpSel, :), 2) / 15;
    
    % Coarsen time discretization
    qInCum = [0; cumsum(qIn)];
    qOutCum = [0; cumsum(qOut)];
    nEl = numel(qInCum);
    iSel = 1:COARSE_SCALE:nEl;
    t = t(iSel(1:end-1));
    dt = (t(2) - t(1));
    qInCum = qInCum(iSel);
    qOutCum = qOutCum(iSel);
    qIn = diff(qInCum);
    qOut = diff(qOutCum);

    switch ACTION
        case DO_CALCULATIONS
            nEl = numel(qOut);
            windowWidth = ceil(INITIAL_WINDOW_WIDTH / STEP);
            pdfWidth = ceil(windowWidth / 1);
            nSteps = ceil((nEl - windowWidth) / STEP);
            muAll = nan(1, nSteps);
            sigmaAll = nan(1, nSteps);
            mcBalanceAll = nan(1, nSteps);
            errAll = nan(1, nSteps);
            qOutConv = zeros(nEl, 1);
            
            % Moisture content balance
            mcBalance = 0;
            
            tic
            for i = 1:STEP:(nEl-windowWidth)
                iAll = ceil(i / STEP);
                mcBalanceAll(iAll) = mcBalance;
                mcBalance = qInCum(i) - qOutCum(i);
                windowWidthLocal = windowWidth;
                err = Inf;
                minErr = Inf;
                while ~(err < ERR_EPSILON) && ~(i + windowWidthLocal == nEl)
                    iSel = i:(i+windowWidthLocal);
                    cumErr = qOut(i) - qOutConv(i);
                    pdfEst = Deconvolution(qIn(iSel), qOut(iSel) - qOutConv(iSel) - cumErr);
                    pdfSel = 1:pdfWidth;
                    tEst = t(i-1+pdfSel)-t(i);
                    [optMu, optSigma, err] = CalcLognormMuSigma(tEst, pdfEst(pdfSel));
                    windowWidthLocal = min(windowWidthLocal * 1.5, nEl - i);
                    if (err < minErr)
                        minErr = err;
                        muAll(iAll) = optMu;
                        sigmaAll(iAll) = optSigma;
                    end
                end
                if (minErr > ERR_EPSILON)
                    if (i == 1)
                        error('Not able to estimate response PDF');
                    else
                        warning('Step %d: didn''t converge', iAll);
                        optMu = nan; muAll(iAll-1);
                        optSigma = nan; sigmaAll(iAll-1);
                        figure(1);
                        plot(t(pdfSel), pdfEst(pdfSel));
                        hold on;
                        plot(t(pdfSel), lognpdfX(t(pdfSel), muAll(iAll), sigmaAll(iAll), dt), 'r');
                        hold off;
                        title(sprintf('mu = %f, sigma = %f, err = %e',  muAll(iAll), sigmaAll(iAll), minErr));
                    end
                end
                iEnd = nEl;
                for j = i:(i+STEP-1)
                    qOutConv(j:iEnd) = qOutConv(j:iEnd) + ...
                        qIn(j) * lognpdfX(t(j:iEnd)-t(j), optMu, optSigma, dt);
                end
                muAll(iAll) = optMu;
                sigmaAll(iAll) = optSigma;
                errAll(iAll) = err; % SumSquares(qOut(i:(i+step-1)), qOutConv(i:(i+step-1)));
                
                % figure(2);
                % hold on;
                % plot(t, qIn/10, 'Color', [0.8, 0.8, 0.8]);
                % plot(t, qOut, 'LineWidth', 2, 'Color', [0.2, 0.2, 0.2]);
                % plot(t, qOutConv, 'LineWidth', 2, 'Color', [0.65, 0.65, 0.65]);
                % plot(t, qOut - qOutConv - cumErr, 'LineStyle', '--', 'Color', [0, 0, 0]);
                % plot([t(i), t(i)], [min(qOut), max(qOut)], 'LineStyle', '-.', 'Color', [0.5, 0.5, 0.5]);
                % hold off;
                % ylim([-1e-5, 3e-4]);
                % xlim([0, t(end)]);
                % legend('Influx (scaled)', 'Original outflux', 'Modelled outflux', ...
                %     'Outflux remaining', 'Modelled until');
                
                fprintf('%f%%\n', iAll / numel(muAll) * 100);
            end
            toc
            
            isNan = isnan(muAll);
            mcBalanceAll(isNan) = [];
            muAll(isNan) = [];
            sigmaAll(isNan) = [];
            errAll(isNan) = [];
            
            save(SAVE_FILE, 'qOutConv', 'mcBalanceAll', 'muAll', 'sigmaAll', 'errAll');
        case LOAD_RESULTS
            load(SAVE_FILE);
    end
    
    close all;
%     iSel = 1:(nSteps * step);
    iSel = ':'; % (errAll <  ERR_EPSILON); % abs(qOut(iSel) - qOutConv(iSel)) < 0.05 * max(qOut(iSel));
%     iSel = sum(reshape(iSel, step, []), 1) > 1;
    fWidth = 390;
    fHeight = 300;
    fH = figure(3);
    set(fH, 'Position', [100, 200, fWidth, fHeight]);
% muAll(mcBalanceAll == max(mcBalanceAll)) = 9.05;
% sigmaAll(mcBalanceAll == max(mcBalanceAll)) = 0.61;
mcBalanceAll = -mcBalanceAll;
i = 2;
while (i <= numel(muAll))
    if ((muAll(i) == muAll(i-1)) && (sigmaAll(i) == sigmaAll(i-1)))
        mcBalanceAll(i) = [];
        muAll(i) = [];
        sigmaAll(i) = [];
        i = i - 1;
    end
    i = i + 1;
end
    [mcBalanceAll, iSort] = sort(mcBalanceAll(iSel)');
    plot(mcBalanceAll, muAll(iSort), '*');
    hold on;
    [muFit, eqStr, muCoeff] = RegressionAnalysis(mcBalanceAll, muAll(iSort), 2);
    plot(mcBalanceAll, muFit, 'r');
    hold off;
    title('Log-normal location parameter (\mu)');
    ylabel('\mu');
    xlabel('Available storage decrement');
    legend('points', eqStr);
    fH = figure(4);
    set(fH, 'Position', [100 + ceil(fWidth * 1.2), 200, fWidth, fHeight]);
    plot(mcBalanceAll, sigmaAll(iSort), '*');
    hold on;
    [sigmaFit, eqStr, sigmaCoeff] = RegressionAnalysis(mcBalanceAll, sigmaAll(iSort), 2);
    plot(mcBalanceAll, sigmaFit, 'r');
    hold off;
    title('Log-normal shape parameter (\sigma)');
    ylabel('\sigma');
    xlabel('Available storage decrement');
    legend('points', eqStr);

    thetaIni = 0;
    vPore = 1;
    
%     % Experiment 1:
%     muEst = 9.795014435936489;
%     sigmaEst = 0.575230541445190;
%     % Functions of log-normal parameters of moisture content
%     fMu = @(x) -1820.71579 * x.^2 + 7.42552 * x + 9.51092;
%     fSigma = @(x) -156.37151 * x.^2 + 5.66496 * x + 0.52888;
    % Experiment 2:
    muEst = 3.32;
    sigmaEst = 0.815;
    fMu = @(x) muCoeff(1) * x.^2 + muCoeff(2) * -x + muCoeff(3);
    fSigma = @(x) sigmaCoeff(1) * x.^2 + sigmaCoeff(2) * -x + sigmaCoeff(3);
    
    LogNormalParams = @(x) {fMu(x), fSigma(x)};
    qOutConvLogn = NumericalConvolution(t, qIn, ...
        @(t, mu, sigma) lognpdfX(t, mu, sigma, dt), muEst, sigmaEst);
    [qOutConvLognVar, dTheta, dMu, dSigma] = NumericalConvolution(t, qIn, ...
        @(t, mu, sigma) lognpdfX(t, mu, sigma, dt), LogNormalParams, thetaIni, vPore);
    
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
    set(axH(2), 'xlim', [t(1), t(end)]);
    set(axH(3), 'xlim', [t(1), t(end)]);
    set(lH{2}(1), 'Color', [0., 0., 0.], 'LineWidth', 1);
    set(lH{2}(2), 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--', 'LineWidth', 2);
    set(lH{2}(3), 'Color', [0.7, 0.7, 0.7], 'LineWidth', 2);
    hgsave(fH, sprintf('./fig/%s_(de)convolution.fig', CASE_NAME));
    
    return

    function pdfEst = Deconvolution(sigIn, sigOut)
        % Fourier transforms of inputs
        sigInF = fft(sigIn);
        sigOutF = fft(sigOut);
        % Deconvolution
        pdfF = sigOutF ./ sigInF;
        % Inverse Fourier transform to get approximation of probability density function
        pdfEst = smooth(ifft(pdfF), 5);
    end
    
    function pdf = lognpdfX(t, mu, sigma, dt)
        pdf = logncdf([t; t(end)+dt], mu, sigma);
        pdf = diff(pdf);
    end
    
    function [optMu, optSigma, optErr] = CalcLognormMuSigma(x, pdf)
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

    function [yFit, eqStr, p_] = RegressionAnalysis(x, y, n)
        if nargin < 3
            n = 2;
        end
        p_ = polyfit(x, y, n);
        yFit = zeros(size(x));
        eqStr = 'y = ';
        for i_ = 1:n+1
            yFit = yFit + p_(i_) * x.^(n+1-i_);
            eqStr = [eqStr, sprintf('%3.2f', p_(i_))];
            if (n+1-i_) == 1
                eqStr = [eqStr, ' * x + '];
            elseif (n+1-i_) > 1
                eqStr = [eqStr, sprintf(' * x^%d + ', n+1-i_)];
            end
        end
%         eqStr(end-1:end) = [];
    end

end