function ModelOutput = ComputeTravelTimes(TimeParams, rainData, rainConcentrationData, ...
        ModelDim, ModelParams)

%% TODO: setting concentration of cloride Comp.consti = XX; Is already concentration
%%       (unit = mol/l). Is very inert.
%%
   
    Const = DefineConstants();
    
    %% Initializing chemical module
    CHEM_MODEL_DIR = '../';
    addpath(genpath(CHEM_MODEL_DIR));
    modus = 0;
    [mIni, Comp, Pm, S, Rp, H] = initialize_ODE(modus, modus);
    ORI = initialize_ORI(Comp, 0);
    
    mInertIni = [1];
    nInertSpecies = numel(mInertIni);
    nReactiveSpecies = numel(mIni);
    mIni = cat(2, mIni, mInertIni);
    nSpecies = numel(mIni);
    iReactiveSpecies = 1:nReactiveSpecies;
    iInertSpecies = nReactiveSpecies:nSpecies;
    iFlushSpecies = [2:5, 8:9, 22];
    nSpeciesDiffuse = numel(iFlushSpecies);
    
    nPhases = 2;
    
    global Call V2 tt Rall
    Call = [];  V2 = []; tt = []; Rall = [];
    options = odeset('OutputFcn','StoreC','Refine',1, 'AbsTol', 1e-8, 'RelTol', 1e-8);
    
    %% Extracting some inputs from input structures
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
    % Volume adjustment factor (Orchestra is calculated for 325 liters of volume), therefore we bring
    % results to volume given in cubic meters in this model.
    pv = zLength .* theta;
    pvAdj = pv / 0.325; % [1, 1];
    
    % Initial concentrations of solutes
    mIni = cat(1, permute(mIni, [3, 1, 2]) * pvAdj(1), ...
        permute(mIni, [3, 1, 2]) * pvAdj(2));
    cIni = mIni ./ repmat(pv, [1, 1, nSpecies]);
    
    % Initialize object to keep information about volumes and concentrations
    % dimensions of arrays (1 x leave time)
    PartInfo = ConcentrationCl(zeros(1, nT), zeros(1, nT, nSpecies));

    %% Initializing output arrays
    qOutTotal = zeros(1, nT);
    mOutTotal = zeros(1, nT, nSpecies);
    cRemaining = nan(nPhases, nT + 1, nSpecies);
    cRemaining(:, 1, :) = cIni;
    mRemaining = nan(nPhases, nT + 1, nSpecies);
    mRemaining(:, 1, :) = mIni;

    %% Main loop
    nPercent = ceil(10 * nT / 100);
    for iT = 1:nT
        tOffset = t(iT);
        tAfter = t(iT:nT) - tOffset;
        
        % Calculate boundaries of log-normally distributed travel times (to save computational
        % resources we calculate only for those times, when something comes out of the drainage
        % system but not until the end of the period).
        tBounds = LogNormalBounds(mu, sigma, Const.NUM_SIGMAS);
        % Select only those time steps that are affected by current injection
        iCalcT = (tAfter <= tBounds(2));
        nCalcT = sum(iCalcT);
        iTend = iT + nCalcT - 1;
        tAfter = tAfter(iCalcT);
        
        if RealGt(rainData(iT), 0, Const.EPSILON)
            % Add (fresh) rainwater to the system. Change volume of liquid and revise mass and
            % concentration of solute in mobile phase
            pvMobUpd = pv(2) + rainData(iT);
            
            %% Recirculation
            DO_RECIRCULATION = false;
            if (DO_RECIRCULATION)
                if ((iT > 1) && (~RealEq(qOutTotal(iT-1), 0, Const.EPSILON)))
                    mRemaining(2, iT, iFlushSpecies) = mRemaining(2, iT, iFlushSpecies) + ...
                        rainData(iT) * mOutTotal(1, iT-1, iFlushSpecies) / qOutTotal(iT-1);
                end
            else
                mRemaining(2, iT, :) = mRemaining(2, iT, :) + ...
                    rainData(iT) * rainConcentrationData(iT);
            end
                
            cRemaining(2, iT, :) = mRemaining(2, iT, :) / pvMobUpd;
            pv(2) = pvMobUpd;
            % Every input impulse of water will cause (log-normal) response at the outlet. All the
            % outflow during a given time step is considered as a particle with unique travel time
            qOutAfter = rainData(iT) * lognpdf(tAfter, mu, sigma) * dt;
            % We integrate volumes of all the particles flowing out at the same time intervals to
            % obtain Leachate volume flux.
            qOutTotal(iT:iTend) = qOutTotal(iT:iTend) + qOutAfter;
            % Increase volume of water leaving the system at future times
            PartInfo = PartInfo.AddVolume(qOutAfter, :, iT:iTend);
        end

        tRange = [tAfter(1), tAfter(1) + dt];

       %% Andre's block
        DO_BIOCHEMISTRY = true;
        if DO_BIOCHEMISTRY
            [~, mChem] = ode15s(@bioreactor, linspace(tRange(1), tRange(2), 3), ...
                mRemaining(1, iT, iReactiveSpecies) / pvAdj(1), options, Comp, Pm, S, Rp, ORI, H);
            mChem = mChem * pvAdj(1);
        else
            mChem = squeeze(repmat(mRemaining(1, iT, iReactiveSpecies), [1, 2, 1]));
        end
        cRemaining(1, iT, iReactiveSpecies) = mChem(end, :) / pv(1);
        mRemaining(1, iT + 1, iReactiveSpecies) = mChem(end, :);
        mRemaining(1, iT + 1, iInertSpecies) = mRemaining(1, iT, iInertSpecies);
        mRemaining(2, iT + 1, :) = mRemaining(2, iT, :);
       %% END
        
        % Solve exchange equation in order to obtain concentrations in both phases at the end of
        % current time step. Different volumes will have different effect on exchange here
        cOutAfter = ConcentrationExchangePhases(tRange, cRemaining(:, iT, iFlushSpecies), kExch, ...
            lambda, pv);
        % Save the remaining mass of solute to output vector
        mRemaining(:, iT + 1, iFlushSpecies) = cOutAfter(:, 2, :) .* repmat(pv, [1, 1, nSpeciesDiffuse]);
        cPartOutR = cat(2, cOutAfter(2, 2, :), PartInfo.GetConcentration(1, iT:iTend, iFlushSpecies));
        pvPartOutR = cat(2, pv(2), PartInfo.GetVolume(1, iT:iTend));
        cPartR = ConcentrationExchangeParticles(...
            [tAfter(1), tAfter(1) + dt], cPartOutR, kExchPart, pvPartOutR);
        PartInfo = PartInfo.SetConcentration(cPartR(2, 2:(nCalcT+1), :), 1, iT:iTend, iFlushSpecies);
        mOutTotal(1, iT, iFlushSpecies) = PartInfo.GetMass(1, iT, iFlushSpecies);
        
        % Withdraw solute leaving together with leachate and compute remaining masses of solutes
        % in both phases
        mRemaining(2, iT + 1, :) = mRemaining(2, iT + 1, :) - mOutTotal(1, iT, :);
        % Remove volume of drained leachate from the volume of liquid in the system.
        pv(2) = pv(2) - qOutTotal(iT);
        % Calculate the remaining concentrations
        cRemaining(:, iT + 1, :) = mRemaining(:, iT + 1, :) ./ repmat(pv, [1, 1, nSpecies]);
        
        % Some checks
        if any(RealLt(reshape(cRemaining(:, iT + 1, iFlushSpecies), 1, []), 0, ...
                Const.CONCENTRATION_EPSILON))
            error('iT = %d: Concentration is negative.', iT);
        end
        
        % Output progress
        if mod(iT, nPercent) == 0
            fprintf('%5.3f%% complete\n', iT / nT * 100);
        end
    end
    
%     fprintf('100%% complete\n\n');
    
    %% Storing outputs to structures
    % Pack outputs to a struct
    ModelOutput = struct();
    ModelOutput.t = t;
    ModelOutput.nT = nT;
    ModelOutput.mIni = mIni;
    ModelOutput.qOutTotal = qOutTotal;
    ModelOutput.mOutTotal = mOutTotal;
    ModelOutput.cRemaining = cRemaining;
    ModelOutput.mRemaining = mRemaining;
    
    ModelOutput.cAll = Call(2:2:end, :);
    
    % Add also concentration at the outlet
    qOutTotalRep = repmat(qOutTotal, [1, 1, nSpecies]);
    ModelOutput.cOutTotal = mOutTotal ./ qOutTotalRep;
    ModelOutput.cOutTotal(qOutTotalRep == 0) = 0;
    
    return
    
    %%
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
                
            if (kExchX == 0)
                cTrend = repmat(cIniX, [1, ntX, 1]);
            else
                if (size(pvX, 3) == 1)
                    pvX = repmat(pvX, [1, 1, nElements]);
                end
                
                vRatioIm = pvX(2, :, :) ./ (pvX(1, :, :) + pvX(2, :, :));
                vRatioM = pvX(1, :, :) ./ (pvX(1, :, :) + pvX(2, :, :));
                
                v01 = vRatioIm .* kExchX;
                v02 = vRatioM .* kExchX;
                v03 = v01 + v02 + lambda;
                v04 = cIniX(2, :, :) .* (v01 - v02 + lambda);
                
                v6 = sqrt((v01 + v02) .^ 2 - lambda .* (2 .* (v01 + v02) + lambda));
                v8 = repmat(-0.5 .* (v03 - v6), [1, ntX, 1]);
                v9 = repmat(-0.5 .* (v03 + v6), [1, ntX, 1]);
                v11 = exp(v8 .* repmat(tX, [1, 1, nElements]));
                v12 = exp(v9 .* repmat(tX, [1, 1, nElements]));
                
                C1_ = 0.5 .* (...
                    cIniX(2, :, :) .* v6 + ...
                    2 .* cIniX(1, :, :) .* v02 + ...
                    v04) ./ v6;
                C2_ = 0.5 .* (...
                    cIniX(2, :, :) .* v6 - ...
                    2 .* cIniX(1, :, :) .* v02 - ...
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
    end

    %%
    function cPart = ConcentrationExchangeParticles(tauX, cIniX, kExchX, pvX)
        % Function to do calculate exchange of solutes between main container and different 
        % particles residing in it. The first element along elements' axis is the main container.
        % All other particles exchange with it, tending towards concentration equilibrium.
        % Input parameters:
        %   tauX        - time span for calculations
        %   cIniX       - initial concentrations of solutes in main container (cIniX(1, 1, :)) and
        %                 other particles (cIniX(2:end, 1, :)). Dimensions of array are:
        %                 nTimeSteps x (nElements + 1) x nSolutes
        %   kExchX      - exchange coefficient between phases
        %   pvX         - volumes of main container (pv(1)) and particles (pv(2:end))
        
        % Get dimensions of a problem
        [~, nElX, nSolutes] = size(cIniX);
        nElX = nElX - 1;
        nTX = numel(tauX);
        
        % Initialize output array
        cPart = zeros(nTX, nElX+1, nSolutes);
        % First element along time axis is initial concentration
        cPart(1, :, :) = cIniX;

        % Find elements with zero volume.
        isPvZero = RealEq(pvX, 0, Const.EPSILON);
        isPvZero(1) = 0;
        isPvNotZero = ~isPvZero;
        isPvNotZero(1) = 0;
        
        if (sum(isPvNotZero) == 0)
            % If all elements have zero volume, they don't exist and we don't need to do 
            % calculations. Just repeat initial concentrations for all time steps
            cPart(2:nTX, 1, :) = repmat(cIniX(1, 1, :), [nTX-1, 1, 1]);
        else
            % Calculate total volume of all elements. This will be needed for computing average
            % concentration
            pvAll = sum(pvX(2:nElX+1));
            % And perform calculations for all time steps
            for iTX = 2:nTX
                tX = tauX(iTX);
                cPart(iTX, 1, :) = cPart(iTX-1, 1, :);
                mElTot = zeros(1, 1, nSolutes);
                dmEl = zeros(1, nElX, nSolutes);
                
                % First we run through all the elements and calculate total masses of solutes in
                % particles and potential mass exchange between them and main container given
                % difference in concentration
                for iEl = 1:nElX
                    if isPvNotZero(iEl+1)
                        mElTot = mElTot + cPart(iTX-1, iEl+1, :) * pvX(iEl+1);
                        dmEl(1, iEl, :) = pvX(1) .* pvX(iEl+1) ./ (pvX(1) + pvX(iEl+1)) .* ...
                            (cPart(iTX-1, 1, :) - cPart(iTX-1, iEl+1, :)) .* ...
                            kExchX * (tX - tauX(iTX-1));
                    end
                end
                % Then calculate mass exchange assuming average concentration of all the elements
                % with the main container
                dmAll = pvX(1) * pvAll / (pvX(1) + pvAll) * ...
                    (cPart(iTX-1, 1, :) - mElTot / pvAll) .* kExchX .* (tX - tauX(iTX-1));
                dmAllSum = sum(dmEl, 2);
                for iSoluteX = 1:nSolutes
                    % For each particle and solute mass exchange is adjusted to match averaged mass 
                    % exchange, but keep relative mass change of particles constant 
                    if ~RealEq(dmAll(1, 1, iSoluteX), 0, Const.EPSILON)
                        if RealEq(dmAllSum(1, 1, iSoluteX), 0, Const.EPSILON)
                            dmAdj = zeros(1, nElX, 1);
                        else
                            dmAdj = dmEl(1, :, iSoluteX) .* dmAll(1, 1, iSoluteX) ./ ...
                                dmAllSum(1, 1, iSoluteX);
                        end
                    else
                        dmAdj = dmEl(1, :, iSoluteX);
                        dmAll(1, 1, iSoluteX) = dmAllSum(1, 1, iSoluteX);
                    end
                    % Update resulting concentrations based on computed before data
                    cPart(iTX, isPvNotZero, iSoluteX) = cPart(iTX-1, isPvNotZero, iSoluteX) + ...
                        dmAdj(isPvNotZero(2:nElX+1)) ./ pvX(isPvNotZero);
                    cPart(iTX, isPvZero, iSoluteX) = 0;
                end
                cPart(iTX, 1, :) = cPart(iTX, 1, :) - dmAll ./ pvX(1);
            end
        end
    end
end