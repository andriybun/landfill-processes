function MainTravelTimes
%% TODO: water must flow even if no rain (two domain?)
%% TODO: mass balance if lambda ~= 0
%%

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
    rainConcentrationData = 0 * ones(size(rainData));
    TimeParams = PrecipitationData.TimeParams;
    
    % Initial concentration (immobile, mobile phases)
    cIni = [1; 1];
    
    %% Some test cases
%     % Case #1:
%     %   initial concentrations = [0; 0]
%     %   constant rain, injection of solute at initial time steps
%     cIni = [0; 0];
%     rainData = 1e-2 * ones(size(rainData));
%     rainConcentrationData(1:5) = 1;
%     % Case #2:
%     %   initial concentrations = [1; 1]
%     %   short clean rain at initial time steps
%     rainData = 0 * rainData;
%     rainData(1:5) = 1e-3;
    %% End test cases
    
    % Time parameters
    tEnd = TimeParams.t(end);   % 100;
    dt = TimeParams.dt;         % 1;
    t = TimeParams.t;           % 0:dt:tEnd;

    % Log-normal parameters
    mu = 0;
    sigma = 1;

    % Other geometrical
    zLength = abs(ModelDim.zin(1) - ModelDim.zin(end));
    % Pore volume
    totalPv = 0.44;
    % Immobile-mobile volume ratio
    beta = 10;
    % Immobile-mobile pore volumes
    theta = [totalPv - totalPv / (1 + beta); totalPv / (1 + beta)];
    pv = zLength .* theta;
    
    % Source/sink rate
    lambda = 0 * 1e-4;
    % Exchange rate between mobile-immobile phases
    kExch = 1e-2;
    % Exchange rate between mobile phase and particles flowing
    kExchPart = -(log(1.5) / exp(mu - sigma^2));
    
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
    
    % Initialize object to keep information about volumes and concentrations
    % dimensions of arrays (entry time x leave time)
    PartInfo = ConcentrationCl(zeros(1, nT), zeros(1, nT));
    
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
        
        if RealGt(rainData(iT), 0, EPSILON)
            % Add (fresh) rainwater to the system. Change volume of liquid and revise mass and
            % concentration of solute in mobile phase
            pvMobUpd = pv(2) + rainData(iT);
            mRemaining(2, iT) = mRemaining(2, iT) + rainData(iT) * rainConcentrationData(iT);
            cRemaining(2, iT) = mRemaining(2, iT) / pvMobUpd;
            pv(2) = pvMobUpd;
            % Every input impulse of water will cause (log-normal) response at the outlet. All the
            % outflow during a given time step is considered as a particle with unique travel time
            qOutAfter = rainData(iT) * lognpdf(tAfter, mu, sigma) * dt;
            % We integrate volumes of all the particles flowing out at the same time intervals to
            % obtain Leachate volume flux.
            qOutTotal(iT:iTend) = qOutTotal(iT:iTend) + qOutAfter;
            % Increase volume of water leaving the system at future times
            PartInfo = PartInfo.AddVolume(qOutAfter, iT:iTend);
        end
        
        % Solve exchange equation in order to obtain concentrations in both phases at the end of
        % current time step. Different volumes will have different effect on exchange here
        cOutAfter = ConcentrationExchangePhases(...
            [tAfter(1), tAfter(1) + dt], cRemaining(:, iT), kExch, lambda, pv);
        % Save the remaining mass of solute to output vector
        mRemaining(:, iT + 1) = cOutAfter(:, 2) .* pv;
        
        cPartOutR = cat(1, repmat(cOutAfter(2, 1), [1, 1, nCalcT]), ...
            reshape(PartInfo.GetConcentration(iT:iTend), 1, 1, []));
        pvPartOutR = cat(1, repmat(pv(2), [1, 1, nCalcT]), ...
            reshape(PartInfo.GetVolume(iT:iTend), 1, 1, []));
        cPartR = ConcentrationExchangePart(...
            [tAfter(1), tAfter(1) + dt], cPartOutR, kExchPart, lambda, pvPartOutR);
        cPart = reshape(cPartR(:, 2, :), [1, nCalcT]);
        PartInfo = PartInfo.SetConcentration(cPart, iT:iTend);
        mOutTotal(iT) = PartInfo.GetMass(iT);
        
        % Withdraw solute leaving together with leachate and compute remaining masses of solutes
        % in both phases
        mRemaining(2, iT + 1) = mRemaining(2, iT + 1) - mOutTotal(iT);
        % Remove volume of drained leachate from the volume of liquid in the system.
        pv(2) = pv(2) - qOutTotal(iT);
        % Calculate the remaining concentrations
        cRemaining(:, iT + 1) = mRemaining(:, iT + 1) ./ pv;
        
        % Some checks
        if any(RealLt(cRemaining(:, iT + 1), 0, EPSILON))
            error('iT = %d: Concentration is negative.', iT);
        end
        if any(RealGt(cRemaining(:, iT + 1), 1, EPSILON))
            error('iT = %d: Concentration is too high.', iT);
        end
        
        % Output progress
        if mod(iT, 1000) == 0
            fprintf('%5.3f%% complete\n', iT / nT * 100);
        end
    end

    toc
    
%     profile off
%     profile viewer
    
    % Results:
    % Out concentration
    cOutRes = mOutTotal(1:nT) ./ qOutTotal(1:nT);
    cOutRes(qOutTotal == 0) = 0;
    % Mass of contaminants removed
    mOutRes = mOutTotal(1:nT);
    % Emission potential
    mRemRes = sum(mRemaining(:, 2:nT+1), 1);

%     plotyy(t, qOutTotal(1:nT), t, mOutTotal(1:nT))
%     return
    
    %% Error check
    if ~RealEq(sum(mIni) - sum(mOutTotal), mRemRes(end), EPSILON)
        warning('ResultCheck:MassBalanceError', 'Absolute error is too high: err = %3.2e', ...
            abs(abs(sum(mIni) - sum(mOutTotal) - mRemRes(end))));
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
        DiffBl.cOutRes = cOutRes(1:nEl) - BaselineRes.cOutRes(1:nEl);
        DiffBl.mOutRes = mOutRes(1:nEl) - BaselineRes.mOutRes(1:nEl);
        DiffBl.mRemRes = mRemRes(1:nEl) - BaselineRes.mRemRes(1:nEl);
        
        fprintf('Error analysis:\n');
        for varIdx = 1:numel(COMP_VARS)
            var = COMP_VARS{varIdx};
            fprintf('\t%s : %f\n', var, max(abs(DiffBl.(var))));
        end
    end
    
    %% Plotting
%     tShow = (TimeParams.daysElapsed > 150) & (TimeParams.daysElapsed < 250);
%     ShowPlots(qOutTotal, mOutTotal, emissionPotential, rainData, lambda, TimeParams, tShow);
    ShowPlots(qOutTotal, mOutTotal, mRemRes, rainData, lambda, TimeParams);
    
    return
    
    function cOut = ConcentrationExchangePhases(tau, cIniX, kExchX, lambda, pv)
        if (tau(1) == tau(end))
            cOut = cIniX;
        else
            cOutAn = Concentration(tau, cIniX, kExchX, lambda, pv);
            cOut = cOutAn;
        end
        
        return
        
        function cTrend = Concentration(tX, cIniX, kExchX, lambda, pvX)
            % Input arrays may have following dimensions:
            %   nPhases x nTimeSteps x nElements (cells, particles etc.)
            
            [nPhases, ~, nElements] = size(cIniX);
            ntX = numel(tX);
            
            vRatioIm = pvX(2, :, :) ./ (pvX(1, :, :) + pvX(2, :, :));
            vRatioM = pvX(1, :, :) ./ (pvX(1, :, :) + pvX(2, :, :));
            
            v01 = vRatioIm .* kExchX;
            v02 = vRatioM .* kExchX;
            v03 = v01 + v02 + lambda;
            v04 = cIniX(2) .* (v01 - v02 + lambda);
            
            v6 = sqrt((v01 + v02) .^ 2 - lambda .* (2 .* (v01 + v02) + lambda));
            v8 = repmat(-0.5 .* (v03 - v6), [1, ntX, 1]);
            v9 = repmat(-0.5 .* (v03 + v6), [1, ntX, 1]);
            v11 = exp(v8 .* repmat(tX, [1, 1, nElements]));
            v12 = exp(v9 .* repmat(tX, [1, 1, nElements]));
            
            C1_ = 0.5 .* (...
                cIniX(2) .* v6 + ...
                2 .* cIniX(1) .* v02 + ...
                v04) ./ v6;
            C2_ = 0.5 .* (...
                cIniX(2) .* v6 - ...
                2 .* cIniX(1) .* v02 - ...
                v04) ./ v6;
            
            C1_ = repmat(C1_, [1, ntX, 1]);
            C2_ = repmat(C2_, [1, ntX, 1]);
            
            cTrend = nan(nPhases, ntX, nElements);
            cTrend(2, :, :) = ...
                C1_ .* v11 + ...
                C2_ .* v12;
            cTrend(1, :, :) = ...
                cTrend(2, :, :) + ...
                (C1_ .* v8 .* v11 + C2_ .* v9 .* v12) ./ repmat(v02, [1, ntX, 1]);
        end
    end

    function cPart = ConcentrationExchangePart(tauX, cIniX, kExchX, lambda, pvX)
        [~, ~, nElements] = size(cIniX);
        ntX = numel(tauX);

        tauX = repmat(tauX, [1, 1, nElements]);
        rate = 1 - exp(kExchX .* tauX);
        cIniX = repmat(cIniX, [1, ntX, 1]);
        cPart = cIniX(2, :, :) + (cIniX(1, :, :) - cIniX(2, :, :)) .* rate;
    end
end