function MainTravelTimesMultiCell
%% TODO: concentration shouldn't drop to zero?
%% TODO: water must flow even if no rain (two domain?)
%% TODO: mass balance if lambda ~= 0
%%

    % Array dimensions:
    %   2d - ([nodes] x [time steps])
    %   3d - ([nodes] x [time steps] x [phases])

    close all

    addpath('../../Common/');
    addpath('../Data/');
    
    % Tolerance
    EPSILON = 1e-10;
    NUM_SIGMAS = 6;
    
    % Dimensions
    zTop = 0;
    zBottom = -1;
    dz = 1;
    ModelDim = InitializeNodes('z', zBottom, zTop, dz);
    ModelDim.zPerc = ModelDim.zin / ModelDim.zin(1);
    
    nZn = ModelDim.znn;
    
    % Load: rainData, TimeParams, StartDate
    PrecipitationData = load('precipitation');
    % Precipitation in meters per time interval dt
    rainData = PrecipitationData.rainData;
    rainConcentrationData = 0 * ones(size(rainData));
    TimeParams = PrecipitationData.TimeParams;
%     rainData = 0 * rainData;
%     rainData(1:2) = 1e-3;
%     rainData(91:92) = 1e-3;
    
    % Time parameters
    tEnd = TimeParams.t(end);   % 100;
    dt = TimeParams.dt;         % 1;
    t = TimeParams.t;           % 0:dt:tEnd;

    % Log-normal parameters
    mu = 0;
    sigma = 1;

    % Pore volume
    theta = [0.4, 0.04];
    pv = permute(ModelDim.dzin * theta, [1, 3, 2]);
    
    % Source/sink rate
    lambda = 0*1e-4;
    % Exchange rate between mobile-immobile phases
    kExch = 1e-2;
    % Exchange rate between mobile phase and particles flowing
    kExchPart = 5e-1;
    
    % Initial concentration (immobile, mobile phases)
    cIni = repmat([1, 1], [nZn, 1]);
    % Initial mass of solute
    mIni = pv .* permute(cIni, [1, 3, 2]);
    
    TimeParams.maxDays = 30;
    nT = TimeParams.maxDays * TimeParams.intervalsPerDay;
    t = t(1:nT);
    
    qOutTotal = zeros(nZn + 1, nT);
    mOutTotal = zeros(nZn, nT);
    cOutTotal = zeros(nZn, nT);
    cRemaining = nan(nZn, nT + 1, 2);
    cRemaining(:, 1, :) = cIni;
    mRemaining = nan(nZn, nT + 1, 2);
    mRemaining(:, 1, :) = mIni;
    
%     profile on
    tic
    
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
        
        % Initialize internodal flux
        qOutAfter = zeros(nZn + 1, nCalcT);
        qOutAfter(nZn + 1, 1) = rainData(iT);
        
        if RealGt(rainData(iT), 0, EPSILON)
            % Every input impulse of water will cause (log-normal) response at the outlet. All the
            % outflow during a given time step is considered as a particle with unique travel time
            muScaled = mu + log(ModelDim.zPerc(1:nZn));
            qOutAfter(1:nZn, :) = rainData(iT) * ...
                lognpdf(repmat(tAfter, [nZn, 1]), repmat(muScaled, [1, nCalcT]), sigma) * dt;
            % We integrate volumes of all the particles flowing out at the same time intervals to
            % obtain Leachate volume flux.
            qOutTotal(:, iT:iTend) = qOutTotal(:, iT:iTend) + qOutAfter(:, :);
        end
        
        % ########### Add comment ###########
        qFluxIn = qOutTotal(:, iT);
        cFluxIn = cat(1, cRemaining(:, iT, 2), rainConcentrationData(iT));
        mFluxIn = qFluxIn .* cFluxIn;
        pvIni = pv;
        pv(:, 1, 2) = pv(:, 1, 2) + qFluxIn(2:nZn+1) - qFluxIn(1:nZn);
        mRemaining(:, iT, 2) = mRemaining(:, iT, 2) + mFluxIn(2:nZn+1) - mFluxIn(1:nZn);
        cRemaining(:, iT, :) = mRemaining(:, iT, :) ./ pv;
        
        % Solve exchange equation in order to obtain concentrations in both phases at the end of
        % current time step. Different volumes will have different effect on exchange here
        cOutAfter = ConcentrationExchangePhases(...
            [tAfter(1), tAfter(1) + dt], cRemaining(:, iT, :), kExch, lambda, pv);
        

%             % We integrate masses of solutes in all particles.
%             % Particles of fresh water also exchange with the surrounding environment and tend to
%             % increase content of solute. The longer particle resides in mobile phase and exchanges 
%             % with it, the higher concentration at outlet will be

        cPart = cat(3, cRemaining(:, iT, 2), cFluxIn(2:nZn+1, 1));
        pvPart = cat(3, pv(:, :, 2), qFluxIn(2:nZn+1, 1));
        cPart = ConcentrationExchangePart(tAfter, cPart, kExchPart, lambda, pvPart);
        mOutTotal(:, iT:iTend) = mOutTotal(:, iT:iTend) + qOutAfter(1:nZn, :) .* cPart;
        
        cOutTotal(:, iT) = mOutTotal(:, iT) ./ qOutTotal(1:nZn, iT);
        isZeroFlux = RealEq(qOutTotal(:, iT), 0, EPSILON);
        cOutTotal(isZeroFlux, iT) = 0;

        % Save the remaining mass of solute to output vector
        mRemaining(:, iT + 1, :) = cOutAfter(:, 2, :) .* pv;
        % Withdraw solute leaving together with leachate and compute remaining masses of solutes
        % in both phases
        mRemaining(:, iT + 1, 2) = mRemaining(:, iT + 1, 2) - mOutTotal(:, iT);
%         [qOutTotal(:, iT), mOutTotal(:, iT), mRemaining(:, iT + 1, 2)]

%         % Add (fresh) rainwater to the system. Change volume of liquid and revise mass and
%         % concentration of solute in the top cell of mobile phase
%         pvMobUpd = pv(:, 1, 2) + qFluxIn(2:nZn+1) - qFluxIn(1:nZn);
%         cRemaining(:, iT, 2) = mRemaining(:, iT, 2) ./ pvMobUpd;
%         pv(:, 1, 2) = pvMobUpd;
        if ~RealEq(sum(sum(pvIni)) + qFluxIn(nZn + 1) - qFluxIn(1), sum(sum(pv)), EPSILON)
            warning('Error on %d-th step\n', iT);
        end

        % Calculate the remaining concentrations
        cRemaining(:, iT + 1, :) = mRemaining(:, iT + 1, :) ./ pv;
        
        % Some checks
        if any(RealLt(cRemaining(:, iT + 1, :), 0, EPSILON))
            error('iT = %d: Concentration is negative.', iT);
        end
        if any(RealGt(cRemaining(:, iT + 1, :), 1, EPSILON))
            error('iT = %d: Concentration is too high.', iT);
        end
        
        % Output progress
        if mod(iT, 1000) == 0
            fprintf('%5.3f%% complete\n', iT / nT * 100);
        end
    end
    
% plot(t, qOutTotal);
% return

    toc
    
%     profile off
%     profile viewer
    
    % Results:
    cOutRes = cOutTotal(1, 1:nT);
    mOutRes = mOutTotal(1, 1:nT);
    mRemRes = sum(sum(mRemaining(:, 2:nT+1, :), 3), 1);

    %% Error check
    if ~RealEq(sum(sum(mIni, 3), 1) - sum(mOutTotal(1, :)), mRemRes(end), EPSILON)
        warning('ResultCheck:MassBalanceError', 'Absolute error is too high: err = %3.2e', ...
            abs(abs(sum(sum(mIni, 3), 1) - sum(mOutTotal(1, :)) - mRemRes(end))));
    end
    %% End error check

    % Validate
    NO_VALIDATION = 0;
    SAVE_RESULTS = 1;
    COMPARE_RESULTS = 2;
    BASELINE_FILE_NAME = '../Data/baseline';
    COMP_VARS = {'cOutRes', 'mOutRes', 'mRemRes'};
    
%     action = SAVE_RESULTS;
    action = COMPARE_RESULTS;
%     action = NO_VALIDATION;
    if (action == SAVE_RESULTS)
        save(BASELINE_FILE_NAME, 'cOutRes', 'mOutRes', 'mRemRes', 'qOutTotal');
    elseif (action == COMPARE_RESULTS)
        BaselineRes = load(BASELINE_FILE_NAME);
        nEl = min(numel(cOutRes), numel(BaselineRes.cOutRes));
% plot(t, cOutRes(1:nEl)); hold on; plot(t, BaselineRes.cOutRes(1:nEl), 'g'); hold off;
% plot(t, qOutTotal(1, 1:nEl)); hold on; plot(t, BaselineRes.qOutTotal(1:nEl), 'g'); hold off;
        DiffBl.cOutRes = cOutRes(1, 1:nEl) - BaselineRes.cOutRes(1:nEl);
        DiffBl.mOutRes = mOutRes(1, 1:nEl) - BaselineRes.mOutRes(1:nEl);
        DiffBl.mRemRes = mRemRes(1, 1:nEl) - BaselineRes.mRemRes(1:nEl);
        
        fprintf('Error analysis:\n');
        for varIdx = 1:numel(COMP_VARS)
            var = COMP_VARS{varIdx};
            fprintf('\t%s : %f\n', var, max(abs(DiffBl.(var))));
        end
    end
    
    %% Plotting
%     tShow = (TimeParams.daysElapsed > 150) & (TimeParams.daysElapsed < 250);
%     ShowPlots(qOutTotal, mOutTotal, emissionPotential, rainData, lambda, TimeParams, tShow);
    ShowPlots(qOutTotal, mOutRes, mRemRes, rainData, lambda, TimeParams);
    
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
    
    function cOut = ConcentrationExchangePhases(tau, cIniX, kExchX, lambda, pv)
        if (tau(1) == tau(end))
            cOut = cIniX;
        else
            cOut = Concentration(tau, cIniX, kExchX, lambda, pv);
        end
        
        return
        
        function cTrend = Concentration(t, cIniX, kExchX, lambda, pv)
            nZnX = size(cIniX, 1);
            nTX = size(t, 2);
            
            vRatioIm = pv(:, :, 2) ./ (pv(:, :, 1) + pv(:, :, 2));
            vRatioM = pv(:, :, 1) ./ (pv(:, :, 1) + pv(:, :, 2));
            
            v01 = vRatioIm .* kExchX;
            v02 = vRatioM .* kExchX;
            v03 = v01 + v02 + lambda;
            v04 = cIniX(:, :, 2) .* (v01 - v02 + lambda);
            
            v6 = sqrt((v01 + v02) .^ 2 - lambda .* (2 .* (v01 + v02) + lambda));
            v8 = repmat(-0.5 .* (v03 - v6), [1, nTX]);
            v9 = repmat(-0.5 .* (v03 + v6), [1, nTX]);
            v11 = exp(v8 .* repmat(t, [nZnX, 1]));
            v12 = exp(v9 .* repmat(t, [nZnX, 1]));
            
            C1_ = 0.5 .* (...
                cIniX(:, :, 2) .* v6 + ...
                2 .* cIniX(:, :, 1) .* v02 + ...
                v04) ./ v6;
            C2_ = 0.5 .* (...
                cIniX(:, :, 2) .* v6 - ...
                2 .* cIniX(:, :, 1) .* v02 - ...
                v04) ./ v6;
            C1_ = repmat(C1_, [1, nTX]);
            C2_ = repmat(C2_, [1, nTX]);

            cTrend = nan(nZnX, nTX, 2);
            cTrend(:, :, 2) = ...
                C1_ .* v11 + ...
                C2_ .* v12;
            cTrend(:, :, 1) = ...
                cTrend(:, :, 2) + ...
                (C1_ .* v8 .* v11 + C2_ .* v9 .* v12) ./ repmat(v02, [1, nTX]);
        end
    end

    function cPart = ConcentrationExchangePart(tau, cIniX, kExchX, lambda, pv)
        cPart = ConcentrationExchangePhases(tau, cIniX, kExchX, 0, pv);
        cPart(:, :, 1) = [];
    end
end