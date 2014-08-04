function ModelOutput = ComputeTravelTimes(TimeParams, RainInfo, ModelDim, ModelParams, prBar)

    DO_BIOCHEMISTRY = ModelParams.DO_BIOCHEMISTRY;
    DO_RECIRCULATION = ModelParams.DO_RECIRCULATION;
    
%% TODO: setting concentration of cloride Comp.consti = XX; Is already concentration
%%       (unit = mol/l). Is very inert.
%%
    
    Const = DefineConstants();
    
    rainData = RainInfo.intensity;
    rainConcentrationData = RainInfo.concentration;
    
	HAS_PROGRESS_BAR = (nargin == 5);
    
    %% Initializing chemical module
    CHEM_MODEL_DIR = '../';
    addpath(genpath(CHEM_MODEL_DIR));
    modus = 0;
    % javaclasspath('/../Orchestra/Interface/OrchestraInterface.jar');
    [mIni, Pm, ORI, ~, ~] = initialize('../MSWS1/Pmatrix/Pmatrix.csv');
    
    mInertIni = ModelParams.mInertIni;
    nInertSpecies = numel(mInertIni);
    nReactiveSpecies = numel(mIni);
    mIni = cat(2, mIni, mInertIni);
    nSpecies = numel(mIni);
    iReactiveSpecies = 1:nReactiveSpecies;
    iInertSpecies = nReactiveSpecies:nSpecies;
    iFlushSpecies = [2:4, 8:9, 24];
    nFlushSpecies = numel(iFlushSpecies);
    
    nPhases = 2;
    
    global Call V2 tt Rall
    Call = [];  V2 = []; tt = []; Rall = [];
    
    %% Extracting some inputs from input structures
	% Time parameters
    tEnd = TimeParams.t(end);
    dt = TimeParams.dt;
    t = TimeParams.t;
    if isfield(TimeParams, 'intervalsPerDay')
        nT = TimeParams.maxDays * TimeParams.intervalsPerDay;
    else
        nT = numel(t);
    end
    t = t(1:nT);
    
    REL_TOL = 1e-5;
    ABS_TOL = 1e-5;
    optionsChem = odeset('InitialStep', dt, ...
        'OutputFcn', ...
        'Store_Orchestra_Results', ...
        'AbsTol', ABS_TOL, ...
        'RelTol', REL_TOL);
    optionsExch = odeset('InitialStep', dt, ...
        'AbsTol', ABS_TOL, ...
        'RelTol', REL_TOL);
    
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
    % pv(2) = 0;
    % Volume of bioreactor model
    vBr = 0.325;
    % Volume adjustment factor 
    pvAdj = pv / vBr;
    
    % Initial concentrations of solutes
    mIni = cat(1, permute(mIni, [3, 1, 2]) * pvAdj(1), ...
        permute(0 * mIni, [3, 1, 2]) * pvAdj(2));
    cIni = mIni ./ repmat(pv, [1, 1, nSpecies]);
    
    % Initialize object to keep information about volumes and concentrations dimensions of arrays 
    % (1 x leave time) + 2 extra particles for water that is staying after tEnd and was originally 
    % in mobile phase
    PartInfo = ConcentrationCl(zeros(1, nT+2), zeros(1, nT+2, nSpecies));
    PartInfo = PartInfo.AddSolute(pv(2), cIni(1, 1, :), 1, nT+2, :);

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
        
        if RealGt(rainData(iT), 0, Const.EPSILON)
            % Add (fresh) rainwater to the system. Change volume of liquid and revise mass and
            % concentration of solute in mobile phase
            pvMobUpd = pv(2) + rainData(iT);
            pv(2) = pvMobUpd;
            % Every input impulse of water will cause (log-normal) response at the outlet. All the
            % outflow during a given time step is considered as a particle with unique travel time
            if iT == nT
                lnPdf = LogNorm.PdfDelayed([tAfter, tAfter+dt]);
                lnPdf = lnPdf(1);
            else
                lnPdf = LogNorm.PdfDelayed(tAfter);
            end
            qOutAfter = rainData(iT) * lnPdf;
            if (DO_RECIRCULATION)
                if ((iT > 1) && (~RealEq(PartInfo.GetVolume(1, iT-1), 0, Const.EPSILON)))
                    cOutPrevStep = mOutTotal(1, iT-1, iFlushSpecies) / PartInfo.GetVolume(1, iT-1);
                    mRemaining(2, iT, iFlushSpecies) = mRemaining(2, iT, iFlushSpecies) + ...
                        rainData(iT) * cOutPrevStep;
                    PartInfo = PartInfo.AddSolute([qOutAfter, rainData(iT) - sum(qOutAfter)], ...
                        repmat(cOutPrevStep, [1, nCalcT+1, 1]), :, [iT:iTend, nT+1], iFlushSpecies);
                end
            else
                % Increase volume of water leaving the system at future times
                PartInfo = PartInfo.AddSolute([qOutAfter, rainData(iT) - sum(qOutAfter)], ...
                    repmat(rainConcentrationData(iT), [1, nCalcT+1, nFlushSpecies]), ...
                    :, [iT:iTend, nT+1], iFlushSpecies);
                mRemaining(2, iT, :) = mRemaining(2, iT, :) + ...
                    rainData(iT) * rainConcentrationData(iT);
            end
            cRemaining(2, iT, :) = mRemaining(2, iT, :) / pvMobUpd;
        end
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
        mRemaining(1, iT + 1, iInertSpecies) = mRemaining(1, iT, iInertSpecies);
        mRemaining(2, iT + 1, :) = mRemaining(2, iT, :);
        % Prepare inputs for exchange equation
        iActivePart = [iT:iTend, nT+1, nT+2];
        cPartOutR = cat(2, cRemaining(1, iT, iFlushSpecies), ...
            PartInfo.GetConcentration(1, iActivePart, iFlushSpecies));
        pvPartOutR = cat(2, pv(1), PartInfo.GetVolume(1, iActivePart));
        % Don't calculate for particles vith (almost) zero volume
        iCalcLog = (pvPartOutR > Const.VOLUME_EPSILON);
        % ... but make sure immobile phase is calculated
        iCalcLog(1) = 1;
        iCalc = find(iCalcLog);
        % Solve exchange equation in order to obtain concentrations in both phases at the end of
        % current time step. Different volumes will have different effect on exchange here
        if (isequal(iCalc, 1))
            cPartR = cPartOutR(:, iCalcLog, :);
        else
            cPartR = MultiConcentrationExchange([tAfter(1), tAfter(1) + dt], ...
                cPartOutR(:, iCalcLog, :), kExch, pvPartOutR(:, iCalcLog));
        end
        % Update computed concentrations for particles
        PartInfo = PartInfo.SetConcentration(cPartR(end, 2:end, :), ...
            1, iT-2+iCalc(2:end), iFlushSpecies);
        % Update masses of solutes per particle
        mOutTotal(1, iT, iFlushSpecies) = PartInfo.GetMass(1, iT, iFlushSpecies);
        % Update total masses of solutes in phases
        mRemaining(1, iT + 1, iFlushSpecies) = cPartR(end, 1, :) * pv(1);
        mRemaining(2, iT + 1, iFlushSpecies) = ...
            sum(PartInfo.GetMass(1, iActivePart(2:end), iFlushSpecies), 2);
        % Remove volume of drained leachate from the volume of liquid in the system.
        pv(2) = pv(2) - PartInfo.GetVolume(1, iT);
        % Calculate the remaining concentrations
        cRemaining(:, iT + 1, :) = mRemaining(:, iT + 1, :) ./ repmat(pv, [1, 1, nSpecies]);
        isPvZero = (pv == 0);
        cRemaining(isPvZero, iT + 1, :) = 0;
        % Some checks
        if any(RealLt(reshape(cRemaining(:, iT + 1, iFlushSpecies), 1, []), 0, ...
                Const.CONCENTRATION_EPSILON))
            error('iT = %d: Concentration is negative.', iT);
        end
        % Output progress
        if (mod(iT, nPercent) == 0) || (iT == nT)
            if HAS_PROGRESS_BAR
                prBar.pvalue = iT / nT * 100;
                %% TODO: redraw GUI (maybe inside set method)
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
    ModelOutput.qOutTotal = PartInfo.GetVolume(1, 1:nT);
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
    
    return
    
    %%
    function cPart = MultiConcentrationExchange(tRange, cIni, kExch, v, Const)
        % Function to calculate exchange of solutes between main container and different particles
        % residing in it. The first element along elements' axis is the main container. This function
        % uses ODE solver.
        % All other particles exchange with it, tending towards concentration equilibrium.
        % Input parameters:
        %   tauX        - time span for calculations
        %   cIniX       - initial concentrations of solutes in main container (cIniX(1, 1, :)) and
        %                 other particles (cIniX(2:end, 1, :)). Dimensions of array are:
        %                 nTimeSteps x (nElements + 1) x nSolutes
        %   kExchX      - exchange coefficient between phases
        %   pvX         - volumes of main container (pv(1)) and particles (pv(2:end))
        
        % Get dimensions of a problem
        [~, nEl, nSolutes] = size(cIni);
        cIni = permute(cIni, [2, 3, 1]);
        cIni = reshape(cIni, [], 1);
        [~, cResRaw] = ode45(@(tX, cX) dC(tX, cX, v', kExch, [1, nEl, nSolutes]), tRange, cIni, ...
            optionsExch);
        cPart = reshape(cResRaw, [], nEl, nSolutes);
        cPart(2:end-1, :, :) = [];
    end
    
    % System of ODE's
    function dCdt = dC(tX, cX, vX, kExch, dimVec)
        nEl = dimVec(2);
        nSolutes = dimVec(3);
        % Resulting array
        dCdt = zeros(nEl * nSolutes, 1);
        % Indices of mobile particles
        iPartMob = 2:nEl;
        % Loop over solutes (this is slightly faster than matrix operations
        for iSolute = 1:nSolutes
            % Concentrations and their changes are stored in a 1D vector, where all solutes are
            % located sequentially. Thus offset of indices for each solute is calculated
            iSoluteOffset = (iSolute - 1) * nEl;
            % Indices of concentrations of current solute in mobile particles
            iPartMobSol = iPartMob + iSoluteOffset;
            % Gradient of concentrations between immobile phase and mobile particles
            gradC = cX(iPartMobSol) - cX(iSoluteOffset + 1);
            % Sum of volumes of immobile phase and mobile particles
            sumVx = vX(1) + vX(iPartMob);
            % Main relationships
            dCdt(iSoluteOffset + 1) = sum(kExch * vX(iPartMob) ./ sumVx .* gradC);
            dCdt(iPartMobSol) = -kExch * vX(1) ./ sumVx .* gradC;
        end
    end
end