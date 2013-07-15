function ModelOutput = ComputeTravelTimes(TimeParams, rainData, rainConcentrationData, ...
        ModelDim, ModelParams, cIni)

    global EPSILON NUM_SIGMAS 
    
	% Time parameters
    tEnd = TimeParams.t(end);
    dt = TimeParams.dt;
    t = TimeParams.t;
    nT = TimeParams.maxDays * TimeParams.intervalsPerDay;
    t = t(1:nT);
    
    % Model params
    mu = ModelParams.mu;
    sigma = ModelParams.sigma;
    totalPv = ModelParams.totalPv;
    beta = ModelParams.beta;
    lambda = ModelParams.lambda;
    kExch = ModelParams.kExch;
    kExchPart = ModelParams.kExchPart;

    % Other geometrical params
    zLength = abs(ModelDim.zin(1) - ModelDim.zin(end));
    theta = [totalPv - totalPv / (1 + beta); totalPv / (1 + beta)];
    pv = zLength .* theta;
    
    % Initial mass of solute
    mIni = cIni .* pv;
    
    % Initialize object to keep information about volumes and concentrations
    % dimensions of arrays (entry time x leave time)
    PartInfo = ConcentrationCl(zeros(1, nT), zeros(1, nT));

    qOutTotal = zeros(1, nT);
    mOutTotal = zeros(1, nT);
    cRemaining = nan(2, nT + 1);
    cRemaining(:, 1) = cIni;
    mRemaining = nan(2, nT + 1);
    mRemaining(:, 1) = mIni;
    
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
    
    fprintf('100%% complete\n\n');
    
    % Pack outputs to a struct
    ModelOutput = struct();
    ModelOutput.t = t;
    ModelOutput.nT = nT;
    ModelOutput.mIni = mIni;
    ModelOutput.qOutTotal = qOutTotal;
    ModelOutput.mOutTotal = mOutTotal;
    ModelOutput.cRemaining = cRemaining;
    ModelOutput.mRemaining = mRemaining;
    
    % Add also concentration at the outlet
    ModelOutput.cOutTotal = mOutTotal ./ qOutTotal;
    ModelOutput.cOutTotal(qOutTotal == 0) = 0;
    
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
        pvX = repmat(pvX, [1, ntX, 1]);
        rate = 1 - exp(kExchX .* tauX);
        cIniX = repmat(cIniX, [1, ntX, 1]);
        cPart = cIniX(2, :, :) + (cIniX(1, :, :) - cIniX(2, :, :)) .* ...
            pvX(1, :, :) ./ (pvX(1, :, :) + pvX(2, :, :)) .* rate;
    end
end