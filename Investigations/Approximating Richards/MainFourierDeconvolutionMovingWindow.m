function MainFourierDeconvolutionMovingWindow
    % This script is designed for analysis of impulse-response systems where the shape of response
    % varies depending on system properties. In particular we investigate influx-outflux system,
    % where the shape of outflux depends on the water storage available insdide the system. This
    % parameter depends on system's history.
    % Inputs:
    %   Parameters that are set within script are described in a section below
    %   SAVE_FILE - *.mat file containing struct with fields: t - time vector, qIn - influx, qOut -
    %               outflux;
    

    %% PARAMETERS
    % Desired tolerance (adjust empirically)
    ERR_EPSILON = 5e-4;
    
    % Time stepping parameters
    COARSE_SCALE = 5;               % Coarsen time stepping by a scale of ...
    STEP = 50;                      % Time step for considering windows
    INITIAL_WINDOW_WIDTH = 16000;   % Number of time steps/intervals per window
    INITIAL_PDF_WIDTH = 10000;
    
    % Ranges of log-normal parameters
    % %% OLD:[7, 10] [0.4, 0.7]
    % %% EX: [0.7, 2.0] [0.9, 1.2]
    % %% ME: [0.7, 2.0] [0.9, 1.2] 
    % %% NL: [2.3, 3.5] [0.8, 1.1]
    % %% UN: [2.3, 3.5] [0.8, 1.1]
    RANGE_MU = [4.0, 7.0];    %               % Range of possible values for log-normal parameter mu
    RANGE_SIGMA = [0.6, 0.9]; %            % Range of possible values for sigma
    DELTA_MU = 2e-2;                        % Step for values of mu within a range
    DELTA_SIGMA = 5e-3;                     % Step for values of sigma within a range
    
    % Do calculations or load saved results
    DO_CALCULATIONS = 1;
    LOAD_RESULTS = 2;
    ACTION = DO_CALCULATIONS;

    %% END PARAMETERS

    close all

    addpath('../../Common');
    
    % Choose for the case to consider
    MAT_FILE_DIR = '../Transforms/mat/';
%     CASE_NAME = 'Richards_Output_Very_Slow_Response.mat';
%     SAVE_FILE = 'DeconvoluteMuSigma_Very_Slow_Response.mat';
%     CASE_NAME = 'ShirishSimulation.mat';
%     SAVE_FILE = 'DeconvoluteMuSigma_ShirishSimulation.mat';
    caseNameVector = {'Extreme', 'Middle_Extreme', 'Netherlands', 'Uniform', ...
        'Repeat_Short_333', 'Repeat_Sin_31', 'Repeat_Medium', 'pulse_zref=-0.50'};
    caseIdx = 4;
    CASE_NAME = sprintf('Richards_Output_rainfall_%s.mat', caseNameVector{caseIdx});
    % Name of savefile
    SAVE_FILE = sprintf('DeconvoluteMuSigma_%s', caseNameVector{caseIdx});
    
    RawData = load([MAT_FILE_DIR CASE_NAME]);
    if ~isfield(RawData, 'qIn') && isfield(RawData, 'fluxIn')
        RawData.qIn = RawData.fluxIn;
    end
    if ~isfield(RawData, 'qOut') && isfield(RawData, 'fluxOut')
        RawData.qOut = RawData.fluxOut;
    end
    if isrow(RawData.t)
        RawData.t = RawData.t';
    end
    if isrow(RawData.qIn)
        RawData.qIn = -RawData.qIn';
    end
    if isrow(RawData.qOut)
        RawData.qOut = -RawData.qOut';
    end
    
    % Get data for processing
    inpSel =  ':';
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

    windowWidth = ceil(INITIAL_WINDOW_WIDTH / STEP);
    pdfWidth = ceil(INITIAL_PDF_WIDTH / STEP);

    %% Determine initial estimate of constant parameters of log-normal distribution
    % iSel = 1:160;
    % CASE_NAME = [CASE_NAME, '_short_fail'];
    % t = t(iSel);
    % qIn = qIn(iSel);
    % qOut = qOut(iSel);
    showApproximation = true;
    pdfEst = Deconvolution(qIn, qOut);
%     pdfEstExp = DeconvolutionExperimental(qIn, qOut);
% plot([pdfEst, pdfEstExp / 30]);
% legend('orig', 'expr');

	[muEst, sigmaEst, ~] = CalcLognormMuSigma(t(1:pdfWidth), pdfEst(1:pdfWidth));
% muEst = 9.4;
% sigmaEst = 0.55;
    fprintf('Initial estimate:\n');
    fprintf('mu = %f\nsigma = %f\n', muEst, sigmaEst);
    confBounds = LogNormalBounds(muEst, sigmaEst, 2);
    fprintf('Recommended PDF width is %d\n', ceil(confBounds(2) / dt) * STEP);
    ComparePdf(t(1:pdfWidth), pdfEst(1:pdfWidth), muEst, sigmaEst, dt);
    qOutConvLogn = NumericalConvolution(t, qIn, ...
        @(t, mu, sigma) lognpdfX(t, mu, sigma, dt), muEst, sigmaEst);
    
    % Plot
    tInfo = struct();
    tInfo.data = t; % / 24 / 60;
    tInfo.axisLabel = 'time [days]';
    
	qInInfo = struct();
    qInInfo.data = qIn;
    qInInfo.name = {'Influx'};
    qInInfo.axisLabel = 'Influx';
    qInInfo.color = [0.75, 0.75, 0.75];
    qInInfo.lineStyle = '-';
    
    qOutInfo = struct();
    if showApproximation
        qOutInfo.data = cat(2, qOut, qOutConvLogn);
        qOutInfo.name = {'Real data', 'Convolution'};
    else
        qOutInfo.data = cat(2, qOut);
        qOutInfo.name = {'Real data'};
    end
    qOutInfo.axisLabel = 'Outflux';
    qOutInfo.color = [0., 0., 0.; ...
                      0.5, 0.5, 0.5];
    qOutInfo.lineStyle = {'-', '--'};
    
    figH = PlotYyCustom(tInfo, qInInfo, qOutInfo, 'NorthEast', ceil(0.9 * [500, 400, 500, 220]));
    hgsave(figH, sprintf('./fig/%s_convolution.fig', CASE_NAME));

    %% Proceed
    
    switch ACTION
        case DO_CALCULATIONS
            nEl = numel(qOut);
            nSteps = ceil((nEl - windowWidth) / STEP);
            muAll = nan(1, nSteps);
            sigmaAll = nan(1, nSteps);
            mcBalanceAll = nan(1, nSteps);
            lastRfAll = zeros(1, nSteps);
            errAll = nan(1, nSteps);
            qOutConv = zeros(nEl, 1);
            
            % Moisture content balance
            mcBalance = 0;
            
            tic
            for i = 1:STEP:(nEl-windowWidth)
                iAll = ceil(i / STEP);
                mcBalanceAll(iAll) = mcBalance;
                mcBalance = qInCum(i) - qOutCum(i);
                lastRfAll(iAll) = i - find(qIn(1:i) > 0, true, 'last');
                windowWidthLocal = windowWidth;
                err = Inf;
                minErr = Inf;
                while ~(err < ERR_EPSILON) && ~(i + windowWidthLocal == nEl)
                    iSel = i:(i+windowWidthLocal);
                    qAdj = qOut(iSel) - qOutConv(iSel);
                    [pks, lcs] = findpeaks(-qAdj);
                    minIdx = find(qAdj == min(qAdj(1:ceil(0.1 * pdfWidth))), true, 'last');
                    if isempty(minIdx)
                        minIdx = 1;
                    end
                    iInc = minIdx;
                    cumErr = -qAdj(minIdx);
                    qAdj = qAdj + cumErr;
                    maxPeak = max(pks(lcs > INITIAL_PDF_WIDTH / STEP));
                    if isempty(maxPeak)
                        minIdx = numel(qAdj);
                    else
                        minIdx = lcs(find(pks == maxPeak, true, 'last'));
                    end
                    pdfEst = Deconvolution(qIn(i+(iInc:minIdx)), qAdj(iInc:minIdx));
                    pdfSel = 1:min(pdfWidth, numel(pdfEst));
                    tEst = t(i-1+pdfSel)-t(i);
                    [optMu, optSigma, err] = CalcLognormMuSigma(tEst, pdfEst(pdfSel));
                    windowWidthLocal = min(windowWidthLocal * 1.5, nEl - i);
%% TODO: select window until the nearest minimum / close to zero
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
                        muAll(iAll) = muAll(iLastNum);
                        sigmaAll(iAll) = sigmaAll(iLastNum);
                        optMu = nan;
                        optSigma = nan;
                        ComparePdf(t(pdfSel), pdfEst(pdfSel), muAll(iAll), sigmaAll(iAll), dt);
                    end
                else
                    iLastNum = iAll;
                end
                iEnd = nEl;
                for j = i:(i+STEP-1)
                    qOutConv(j:iEnd) = qOutConv(j:iEnd) + ...
                        qIn(j) * lognpdfX(t(j:iEnd)-t(j), muAll(iAll), sigmaAll(iAll), dt);
                end
                muAll(iAll) = optMu;
                sigmaAll(iAll) = optSigma;
                errAll(iAll) = err;
                
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
            lastRfAll(isNan) = [];
            muAll(isNan) = [];
            sigmaAll(isNan) = [];
            errAll(isNan) = [];
            
            save(SAVE_FILE, 'qOutConv', 'mcBalanceAll', 'muAll', 'sigmaAll', 'errAll');
        case LOAD_RESULTS
            load(SAVE_FILE);
    end
    
    %% Plotting
    iSel = ':';
    mcBalanceAll = mcBalanceAll;
% mcBalanceAll = cat(2, mcBalanceAll, 0.018);
% muAll = cat(2, muAll, 9.05);
% sigmaAll = cat(2, sigmaAll, 0.61);
    [mcBalanceAll, iSort] = sort(mcBalanceAll(iSel)');


    close all;
    fWidth = 390;
    fHeight = 300;
    fH = figure(3);
    set(fH, 'Position', [100, 200, fWidth, fHeight]);
    plot(mcBalanceAll, muAll(iSort), '*');
    hold on;
    [muFit, eqStr, muCoeff] = RegressionAnalysis(mcBalanceAll, muAll(iSort), 2);
%muCoeff = [-2539.64, -20.79, 9.51];
    plot(mcBalanceAll, muFit, 'r');
    hold off;
    title('Log-normal location parameter (\mu)');
    ylabel('\mu');
    xlabel('Decrease in available storage');
    legend('points', eqStr);
    fH = figure(4);
    set(fH, 'Position', [100 + ceil(fWidth * 1.2), 200, fWidth, fHeight]);
    plot(mcBalanceAll, sigmaAll(iSort), '*');
    hold on;
    [sigmaFit, eqStr, sigmaCoeff] = RegressionAnalysis(mcBalanceAll, sigmaAll(iSort), 2);
%sigmaCoeff = [-95.53, -4.42, 0.55];
    plot(mcBalanceAll, sigmaFit, 'r');
    hold off;
    title('Log-normal shape parameter (\sigma)');
    ylabel('\sigma');
    xlabel('Decrease in available storage');
    legend('points', eqStr);

    thetaIni = 0;
    vPore = 1;
    
    
    
    
    % Minimum storage capacity change
    fMu = @(x) muCoeff(1) * max(-x, mcBalanceAll(1)).^2 + ...
        muCoeff(2) * max(-x, mcBalanceAll(1)) + muCoeff(3);
    fSigma = @(x) sigmaCoeff(1) * max(-x, mcBalanceAll(1)).^2 + ...
        sigmaCoeff(2) * max(-x, mcBalanceAll(1)) + sigmaCoeff(3);
    
    LogNormalParams = @(x) {fMu(x), fSigma(x)};
    qOutConvLogn = NumericalConvolution(t, qIn, ...
        @(t, mu, sigma) lognpdfX(t, mu, sigma, dt), muEst, sigmaEst);
    [qOutConvLognVar, dTheta, dMu, dSigma] = NumericalConvolution(t, qIn, ...
        @(t, mu, sigma) lognpdfX(t, mu, sigma, dt), LogNormalParams, thetaIni, vPore);
    
    % Scale time
    t = t / 24 / 60;
    
    % Plot
    fH = figure(2);
    set(fH, 'Position', [500, 400, 700, 350]);
    plotyy(t, qIn, t, cat(2, qOutConvLogn, qOutConvLognVar, qOut));
    legend('Influx', 'Convolution (lognormal)', 'Convolution (lognormal, estimated dependency)', ...
        'Real data', 'Location', 'Northwest');
    xlabel('Time, days');
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
        EPS = 1e-14;
        isZero = RealEq(real(sigInF), 0, EPS) & RealEq(imag(sigInF), 0, EPS);
        pdfF(isZero) = complex(0.0, 0.0);
        % Inverse Fourier transform to get approximation of probability density function
        pdfEst = smooth(ifft(pdfF), 5);
    end

%% Experimental

    function pdfEst = DeconvolutionExperimental(sigIn, sigOut)
        % Fourier transforms of inputs
        sigInF = fft(AddExtraZeros(diff(sigIn)));
        sigOutF = fft(AddExtraZeros(diff(sigOut)));
        % Deconvolution
        pdfF = sigOutF ./ sigInF;
        EPS = 1e-14;
        isZero = RealEq(real(sigInF), 0, EPS) & RealEq(imag(sigInF), 0, EPS);
        pdfF(isZero) = complex(0.0, 0.0);
        % Inverse Fourier transform to get approximation of probability density function
        pdfEst = smooth(cumsum(RemoveExtraZeros(ifft(pdfF))), 5);
    end

    function dOut = AddExtraZeros(dIn)
        dOut = cat(1, 0, dIn, [0; 0; 0]);
    end

    function dOut = RemoveExtraZeros(dIn)
        dOut = dIn(2:(end-2));
    end
%% END Experimental

    function pdf = lognpdfX(t, mu, sigma, dt)
        pdf = logncdf([t; t(end)+dt], mu, sigma);
        pdf = diff(pdf);
    end

    function [optMu, optSigma, optErr] = CalcLognormMuSigma(x, pdf)
       if any(isnan(pdf))
            optMu = nan;
            optSigma = nan;
            optErr = nan;
            return
        end
        pdf = max(0, pdf);
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

    function ComparePdf(t_, pdfEst_, mu_, sigma_, dt_)
        figH = figure(1);
        set(figH, 'position', ceil(0.9 * [400, 300, 350, 220]));
        plot(t_, pdfEst_, 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--');
        hold on;
        pdfAnalytic = lognpdfX(t_, mu_, sigma_, dt_);
        if showApproximation
            plot(t_, pdfAnalytic, 'Color', [0.2, 0.2, 0.2]);
        end
        maxVal = max(max(pdfAnalytic), max(pdfEst));
        ylim([-0.1 * maxVal, maxVal]);
        hold off;
        if showApproximation
            legend('Estimated PDF', 'Log-normal');
            err_ = SumSquares(max(0, pdfEst_), pdfAnalytic);
            title(sprintf('%s = %2.3f, %s = %2.3f, err = %2.2e',  ...
                '\mu', mu_, '\sigma', sigma_, err_));
        else
            legend('Estimated PDF');
        end
        hgsave(figH, sprintf('./fig/%s_PDF_est.fig', CASE_NAME));
    end
    
end