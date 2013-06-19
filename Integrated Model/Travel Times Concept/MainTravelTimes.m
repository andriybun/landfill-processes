function MainTravelTimes
    close all

    addpath('../../Common/');
    addpath('../Data/');
    
    % Tolerance
    EPSILON = 1e-10;
    NUM_SIGMAS = 6;
    
    % Dimensions
    zTop = 0;
    zBottom = -1;
    dz = 0.05;
    ModelDim = InitializeNodes('z', zBottom, zTop, dz);
    
    % Load: rainData, TimeParams, StartDate
    PrecipitationData = load('precipitation');
    % Precipitation in meters per time interval dt
    rainData = PrecipitationData.rainData;
    TimeParams = PrecipitationData.TimeParams;
%     rainData = 0 * rainData;
%     rainData(1:2) = 1e-3;
%     rainData(91:92) = 1e-3;
    
%     rainData = repmat(rainData, [1, 2]);
%     TimeParams.maxDays = TimeParams.maxDays * 2;
%     TimeParams.numIntervals = TimeParams.maxDays * TimeParams.intervalsPerDay;
%     TimeParams.t = 0:TimeParams.dt:TimeParams.maxDays - EPSILON;
%     TimeParams.daysElapsed = TimeParams.t;
    
    % Time parameters
    tEnd = TimeParams.t(end);   % 100;
    dt = TimeParams.dt;         % 1;
    t = TimeParams.t;           % 0:dt:tEnd;

    % Log-normal parameters
    mu = 0;
    sigma = 0.999;

    % Other geometrical
    zLength = abs(ModelDim.zin(1) - ModelDim.zin(end));
    % Pore volume
    theta = [0.4; 0.04];
    pv = zLength .* theta;
    
    % Source/sink rate
    lambda = 0 * 1e-4;
    % Exchange rate between mobile-immobile phases
    kExch = 1e-2;
    % Diffusion constant for particles
    diffConst = 7e-1;
    
    % Initial concentration (immobile, mobile phases)
    cIni = [1; 0];
    % Initial mass of solute
    mIni = cIni .* pv;
    
%     TimeParams.maxDays = 30;
    nT = TimeParams.maxDays * TimeParams.intervalsPerDay;
    t = t(1:nT);
    
    qOutTotal = zeros(1, nT);
    mOutTotal = zeros(1, nT);
    cRemaining = nan(2, nT + 1);
    cRemaining(:, 1) = cIni;
    mRemaining = nan(2, nT + 1);
    mRemaining(:, 1) = mIni;
    
%     profile on
    tic
    
    % mRemaining = mIni;
    for iT = 1:nT
        tOffset = t(iT);
        tAfter = t(iT:nT) - tOffset;
        
        % Calculate boundaries of log-normally distributed travel times (to save computational
        % resources we calculate only for those times, when something comes out of the drainage
        % system but not until the end of the period).
        tBounds = LogNormalBounds(mu, sigma, NUM_SIGMAS);
        % Select only those time steps that are affected by current injection
        iCalcT = (tAfter <= tBounds(2));
        nCalcT = sum(iCalcT);
        iTend = iT + nCalcT - 1;
        tAfter = tAfter(iCalcT);
        
        if RealGt(rainData(iT), 0, EPSILON)
            % Add fresh rainwater to the system. Change volume of liquid and revise concentration
            % in mobile phase
            pvMobUpd = pv(2) + rainData(iT);
            cRemaining(2, iT) = cRemaining(2, iT) * pv(2) / pvMobUpd;
            pv(2) = pvMobUpd;
            
            % Every input impulse of water will cause (log-normal) response at the outlet. All the
            % outflow during a given time step is considered as a particle with unique travel time.
            qOutAfter = rainData(iT) * lognpdf(tAfter, mu, sigma) * dt;
            % We integrate volumes of all the particles flowing out at the same time intervals to
            % obtain Leachate volume flux.
            qOutTotal(iT:iTend) = qOutTotal(iT:iTend) + qOutAfter;
        end
        
%         % Initial masses of solute(s) available in the system
%         mStepIni = cRemaining(:, iT) .* pv;
%         % Exchange equation is solved with respect to masses of solutes, since volumes of different
%         % phases differ. Thus different volumes will have different effect on exchange
        vRatio = pv(1) / pv(2);
        cOutAfter = OutConcentration([tAfter(1), tAfter(1) + dt], cRemaining(:, iT), kExch, lambda, vRatio);
%         plot(cOutAfter');
%         mOutAfter = OutConcentration([tAfter(1), tAfter(1) + dt], mStepIni, kExch, lambda);
%         % Convert obtained masses of solutes after exchange to concentrations
%         cOutAfter = mOutAfter ./ [pv, pv];
%         % Save the remaining mass of solute to output vector
%         mRemaining(:, iT + 1) = mOutAfter(:, 2);
        mRemaining(:, iT + 1) = cOutAfter(:, 2) .* pv;
        
        if RealGt(rainData(iT), 0, EPSILON)
            % We integrate masses of solutes in all particles.
            % Particles of fresh water also exchange with the surrounding environment and tend to
            % increase content of solute. The longer particle resides in mobile phase and exchanges 
            % with it, the higher concentration at outlet will be
%             cPart = cOutAfter(2, 2) * (1 - exp(-diffConst * tAfter));
            a = 3e+0;
            b = 3e+0;
            cPart = cOutAfter(2, 2) * exp(-a ./ (tAfter .^ b));
            mOutTotal(iT:iTend) = mOutTotal(iT:iTend) + qOutAfter .* cPart;
        end
        
        % Withdraw solute leaving together with leachate and compute remaining masses of solutes
        % in both phases
        mRemaining(2, iT + 1) = mRemaining(2, iT + 1) - mOutTotal(iT);
        % Remove volume of drained leachate from the volume of liquid in the system.
        pv(2) = pv(2) - qOutTotal(iT);
        % Calculate the remaining concentrations
        cRemaining(:, iT + 1) = mRemaining(:, iT + 1) ./ pv;
        
%         % Some checks
%         if any(cRemaining(:, iT + 1) < 0)
%             error('iT = %d: Concentration is negative.', iT);
%         end
%         if any(cRemaining(:, iT + 1) > 1)
%             error('iT = %d: Concentration is too high.', iT);
%         end
        
        % Output progress
        if mod(iT, 1000) == 0
            fprintf('%5.3f%% complete\n', iT / nT * 100);
        end
    end

%     plot(qOutTotal, mOutTotal(1:nT) ./ qOutTotal(1:nT), 'x');
%     return
    
    toc
    
%     profile off
%     profile viewer
    
    %% Error check
    mEnd = cRemaining(:, end) .* pv;
    if ~RealEq(sum(mIni) - sum(mOutTotal), sum(mEnd), EPSILON)
        warning('ResultCheck:MassBalanceError', 'Absolute error is too high: err = %3.2e', ...
            abs(abs(sum(mIni) - sum(mOutTotal) - sum(mEnd))));
    end
    %% End error check
    
    % Validate
    NO_VALIDATION = 0;
    SAVE_RESULTS = 1;
    COMPARE_RESULTS = 2;
    BASELINE_FILE_NAME = '../Data/baseline';
    COMP_VARS = {'cOutRes', 'mOutRes', 'cRemRes'};
    
    % Results:
    % Out concentration
    cOutRes = mOutTotal(1:nT) ./ qOutTotal(1:nT);
    cOutRes(qOutTotal == 0) = 0;
    mOutRes = mOutTotal(1:nT);
    cRemRes = sum(cRemaining(:, 1:nT));
    emissionPotential = sum(mRemaining(:, 1:nT), 1);
    
%     action = SAVE_RESULTS;
    action = COMPARE_RESULTS;
    if (action == SAVE_RESULTS)
        save(BASELINE_FILE_NAME, 'cOutRes', 'mOutRes', 'cRemRes', 'qOutTotal');
    elseif (action == COMPARE_RESULTS)
        BaselineRes = load(BASELINE_FILE_NAME);
        nEl = min(numel(cOutRes), numel(BaselineRes.cOutRes));
        DiffBl.cOutRes = cOutRes(1:nEl) - BaselineRes.cOutRes(1:nEl);
        DiffBl.mOutRes = mOutRes(1:nEl) - BaselineRes.mOutRes(1:nEl);
        DiffBl.cRemRes = cRemRes(1:nEl) - BaselineRes.cRemRes(1:nEl);
        
        fprintf('Error analysis:\n');
        for varIdx = 1:numel(COMP_VARS)
            var = COMP_VARS{varIdx};
            fprintf('\t%s : %f\n', var, max(abs(DiffBl.(var))));
        end
    end
    
    %% Plotting
%     tShow = (TimeParams.daysElapsed > 150) & (TimeParams.daysElapsed < 250);
%     ShowPlots(qOutTotal, mOutTotal, emissionPotential, rainData, lambda, TimeParams, tShow);
    ShowPlots(qOutTotal, mOutTotal, emissionPotential, rainData, lambda, TimeParams);
    
%     %% Compare with fourier transform
%     lognpdfVec = lognpdf(t, mu, sigma) * dt;
%     
%     rainDataF = fft(rainData);
%     lognpdfVecF = fft(lognpdfVec);
%     
%     qOutF = rainDataF .* lognpdfVecF;
%     qOutIft = ifft(qOutF);
% 
%     plot(t, cat(1, qOutTotal, qOutIft));
%     legend('Numerical integration', 'Fourier transform');

    return
    
    function cOut = OutConcentration(tau, cRemaining, kExch, lambda, vRatio)
        if (tau(1) == tau(end))
            cOut = cRemaining;
        else
            cOutAn = Concentration(tau, cRemaining, kExch, lambda, vRatio);
            cOut = cOutAn;
        end
        
        return
        
        function cTrend = Concentration(t, cIni, kExch, lambda, vRatio)
            v1 = vRatio .^ 4 .* kExch .^ 2;
            v2 = 2 .* vRatio .^ 2 .* kExch .^ 2;
            v3 = 2 .* vRatio .^ 3 .* kExch .* lambda;
            v4 = 2 .* vRatio .* kExch .* lambda;
            v5 = lambda .^ 2 .* vRatio .^ 2;
            
            C1_ = 0.5 .* (cIni(2) .* ...
                sqrt(v1 + ...
                v2 - ...
                v3 + ...
                kExch .^ 2 + ...
                v4 + ...
                v5) + ...
                2 .* cIni(1) .* vRatio .^ 2 .* kExch - ...
                cIni(2) .* vRatio .^ 2 .* kExch + ...
                cIni(2) .* kExch + ...
                cIni(2) .* lambda .* vRatio) ./ ...
                sqrt(v1 + ...
                v2 - ...
                v3 + ...
                kExch .^ 2 + ...
                v4 + ...
                v5);
            C2_ = -0.5 .* (2 .* cIni(1) .* vRatio .^ 2 .* kExch - ...
                cIni(2) .* vRatio .^ 2 .* kExch + ...
                cIni(2) .* kExch + ...
                cIni(2) .* lambda .* vRatio - ...
                cIni(2) .* ...
                sqrt(v1 + ...
                v2 - ...
                v3 + ...
                kExch .^ 2 + ...
                v4 + ...
                v5)) ./ ...
                sqrt(v1 + ...
                v2 - ...
                v3 + ...
                kExch .^ 2 + ...
                v4 + ...
                v5);
            
            cTrend = nan(2, numel(t));
            cTrend(1, :) = (0.5 .* C2_ - 0.5 .* C2_ ./ vRatio .^ 2 + ...
                (-0.5 .* C2_ .* lambda ./ vRatio - ...
                0.5 .* C2_ .* ...
                sqrt(v1 + ...
                v2 - ...
                v3 + ...
                kExch .^ 2 + ...
                v4 + ...
                v5) ./ vRatio .^ 2) ./ kExch) .* ...
                exp(-0.5 .* (vRatio .^ 2 .* kExch + kExch + lambda .* vRatio + ...
                sqrt(v1 + ...
                v2 - ...
                v3 + ...
                kExch .^ 2 + ...
                v4 + ...
                v5)) .* t ./ vRatio) + ...
                (0.5 .* C1_ - 0.5 .* C1_ ./ vRatio .^ 2 + ...
                (-0.5 .* C1_ .* lambda ./ vRatio + ...
                0.5 .* C1_ .* sqrt(v1 + ...
                v2 - ...
                v3 + ...
                kExch .^ 2 + ...
                v4 + ...
                v5) ./ vRatio .^ 2) ./ kExch) .* exp(0.5 .* ...
                (-vRatio .^ 2 .* kExch - ...
                kExch-lambda .* vRatio + ...
                sqrt(v1 + ...
                v2 - ...
                v3 + ...
                kExch .^ 2 + ...
                v4 + ...
                v5)) .* t ./ vRatio);
            cTrend(2, :) = C1_ .* exp(0.5 .* (-vRatio .^ 2 .* kExch - ...
                kExch-lambda .* vRatio + ...
                sqrt(v1 + ...
                v2 - ...
                v3 + ...
                kExch .^ 2 + ...
                v4 + ...
                v5)) .* t ./ vRatio) + ...
                C2_ .* exp(-0.5 .* (vRatio .^ 2 .* kExch + kExch + ...
                lambda .* vRatio + ...
                sqrt(v1 + ...
                v2 - ...
                v3 + ...
                kExch .^ 2 + ...
                v4 + ...
                v5)) .* t ./ vRatio);
%             v1 = lambda^2 + 4 * kExch^2;
%             v2 = sqrt(v1);
%             v3 = lambda + 2 * kExch;
%             v4 = exp(0.5 * (-v3 + v2) * t);
%             v5 = exp(-0.5 * (v3 + v2) * t);
%             
%             cv1 = cIni(2) * lambda;
%             cv2 = 2 * cIni(1) * kExch;
%             cv3 = cIni(2) * sqrt(v1);
%             
%             C1_ = 0.5 * (cv1 + cv2 + cv3) / v2;
%             C2_ = 0.5 * (-cv1 - cv2 + cv3) / v2;
%             
%             cTrend = nan(2, numel(t));
%             cTrend(1, :) = -0.5 / kExch * (...
%                 C1_ * v4 * lambda - C1_ * v4 * v2 + ...
%                 C2_ * v5 * lambda + C2_ * v5 * v2);
%             cTrend(2, :) = C1_ * v4 + C2_ * v5;
        end
    end
end