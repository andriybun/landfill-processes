function ModelOutput = Compute(TimeParams, RainInfo, ModelDim, ModelParams, prBar)

    DO_BIOCHEMISTRY = ModelParams.DO_BIOCHEMISTRY;
    Const = DefineConstants();
    rainData = RainInfo.intensity;
    
	HAS_PROGRESS_BAR = (nargin == 5);
    if HAS_PROGRESS_BAR
        prBar.pvalue = 0;
    end
    
    %% Initializing chemical module
    CHEM_MODEL_DIR = '../';
    addpath(genpath(CHEM_MODEL_DIR));
    [mIni, Pm, ORI, ~, ~] = initialize('../MSWS1/Pmatrix/Pmatrix.csv');
    
    mInertIni = ModelParams.mInertIni;
    nInertSpecies = numel(mInertIni);
    nReactiveSpecies = numel(mIni);
    % Initial masses of compounds (1 x [num species])
    mIni = cat(2, mIni, mInertIni);
    nSpecies = numel(mIni);
    iReactiveSpecies = 1:nReactiveSpecies;
    iInertSpecies = nReactiveSpecies:nSpecies;
    iFlushSpecies = [2:4, 8:9, 24];
    nFlushSpecies = numel(iFlushSpecies);
    
    SpeciesInfo = struct(...
        'n', nSpecies, ...
        'iInert', iInertSpecies, ...
        'iFlush', iFlushSpecies, ...
        'nFlush', nFlushSpecies);
    
    nPhases = 2;
    
    global Call V2 tt Rall
    Call = [];  V2 = []; tt = []; Rall = [];
    
    %% Extracting some inputs from input structures
	% Time parameters
    tEnd = TimeParams.t(end);
    dt = TimeParams.dt;
    t = TimeParams.t;
	nT = numel(t);
    
    Const.REL_TOL = 1e-5;
    Const.ABS_TOL = 1e-5;
    optionsChem = odeset('InitialStep', dt, ...
        'OutputFcn', ...
        'Store_Orchestra_Results', ...
        'AbsTol', Const.ABS_TOL, ...
        'RelTol', Const.REL_TOL);
    
    % Model params
    LogNorm = ModelParams.LogNorm;
    totalPv = ModelParams.totalPv;
    beta = ModelParams.beta;
    lambda = ModelParams.lambda;
    kExch = ModelParams.kExch;

    % Other geometrical params
    zLength = abs(ModelDim.zin(1) - ModelDim.zin(end));
    theta = [totalPv - totalPv / (1 + beta); totalPv / (1 + beta)];
    % Volume adjustment factor (Orchestra is calculated for 325 liters of volume), therefore we bring
    % results to volume given in cubic meters in this model.
    pv = zLength .* theta;
    pv(2) = 1e-3;
    % pv(2) = 0;
    % Volume of bioreactor model
    vBr = 0.325;
    % Volume adjustment factor 
    pvAdj = pv / vBr;
    
    % Initial mass - add extra dimension for mobile phase and permute:
    % [num phases] x 1 x [num species]
    mIni = cat(1, permute(mIni, [3, 1, 2]) * pvAdj(1), ...
        permute(mIni, [3, 1, 2]) * pvAdj(2));
    % Initial concentrations of solutes
    cIni = mIni ./ repmat(pv, [1, 1, nSpecies]);

    % Some indices for the following particles
    nPartGeneral = 2;
    iImmobile = 1;
    iMobileIni = 2;
    % Initialize object to keep information about volumes and concentrations dimensions of arrays 
    % (1 x leave time) + 2 extra particles for water that is staying after tEnd and was originally 
    % in mobile phase (located in the beginning of array)
    PartInfo = ConcentrationCl(zeros(1, nT + nPartGeneral), zeros(1, nT + nPartGeneral, nSpecies));
    % Immobile phase
    PartInfo = PartInfo.AddSolute(pv(1), cIni(1, 1, :), 1, iImmobile, :);
    % Originally present in mobile phase
    PartInfo = PartInfo.AddSolute(pv(2), cIni(2, 1, :), 1, iMobileIni, :);

    %% Initializing output arrays
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
        tBounds = LogNorm.bounds;
        % Select only those time steps that are affected by current injection
        iCalcT = (tAfter <= tBounds(2));
        nCalcT = sum(iCalcT);
        iTend = iT + nCalcT - 1;
        tAfter = tAfter(iCalcT);
        
        % Apply rainfall
        [pv, PartInfo, mRemaining] = Rainfall(...
            t, dt, iT, pv, mRemaining, RainInfo, PartInfo, SpeciesInfo, ModelParams, Const);
        if RealGt(rainData(iT), 0, Const.EPSILON)
            cRemaining(2, iT, :) = mRemaining(2, iT, :) /  pv(2);
        end
        
        % Apply biochemistry
        if DO_BIOCHEMISTRY
            tRange = [tAfter(1), tAfter(1) + dt];
            [~, mChem] = ode15s(@bioreactor, tRange, ...
                mRemaining(1, iT, iReactiveSpecies) / pvAdj(1), optionsChem, Pm, ORI);
            mChem = mChem * pvAdj(1);
        else
            mChem = squeeze(repmat(mRemaining(1, iT, iReactiveSpecies), [1, 2, 1]));
        end
        cRemaining(1, iT, iReactiveSpecies) = mChem(end, :) / pv(1);
        mRemaining(1, iT + 1, iReactiveSpecies) = mChem(end, :);
        
        [pv, PartInfo, mRemaining, cRemaining] = ExchangePhases(...
            t, dt, iT, pv, mRemaining, cRemaining, PartInfo, SpeciesInfo, ModelParams, Const);
        
        % Remove volume of drained leachate from the volume of liquid in the system.
        pv(2) = pv(2) - PartInfo.GetVolume(1, nPartGeneral + iT);
        % Calculate the remaining concentrations
        cRemaining(:, iT + 1, :) = mRemaining(:, iT + 1, :) ./ repmat(pv, [1, 1, SpeciesInfo.n]);
        isPvZero = (pv == 0);
        cRemaining(isPvZero, iT + 1, :) = 0;
        % Update masses of solutes per particle
        mOutTotal(1, iT, SpeciesInfo.iFlush) = PartInfo.GetMass(1, nPartGeneral + iT, ...
            SpeciesInfo.iFlush);
        
        % Some checks
        if any(RealLt(reshape(cRemaining(:, iT + 1, iFlushSpecies), 1, []), 0, ...
                Const.CONCENTRATION_EPSILON))
            error('iT = %d: Concentration is negative.', iT);
        end
        % Output progress
        if (mod(iT, nPercent) == 0) || (iT == nT)
            if HAS_PROGRESS_BAR
                prBar.pvalue = iT / nT * 100;
            else
                fprintf('%5.3f%% complete\n', iT / nT * 100);
            end
        end
    end
    
    %% Storing outputs to structures
    % Pack outputs to a struct
    ModelOutput = struct();
    ModelOutput.t = t;
    ModelOutput.nT = nT;
    ModelOutput.mIni = mIni;
    ModelOutput.qIn = rainData;
    ModelOutput.qOutTotal = PartInfo.GetVolume(1, nPartGeneral + (1:nT));
    ModelOutput.mOutTotal = mOutTotal;
    ModelOutput.cRemaining = cRemaining;
    ModelOutput.mRemaining = mRemaining;
    ModelOutput.cAll = Call(2:2:end, :);
    
    % Add also concentration at the outlet
    qOutTotalRep = repmat(ModelOutput.qOutTotal, [1, 1, nSpecies]);
    ModelOutput.cOutTotal = mOutTotal ./ qOutTotalRep;
    isFluxZero = RealEq(ModelOutput.qOutTotal, 0, Const.VOLUME_EPSILON);
    ModelOutput.cOutTotal(:, isFluxZero, :) = ...
        (pv(1) * ModelOutput.cRemaining(1, isFluxZero, :) + ...
        pv(2) * ModelOutput.cRemaining(2, isFluxZero, :)) ./ (pv(1) + pv(2));
    
end