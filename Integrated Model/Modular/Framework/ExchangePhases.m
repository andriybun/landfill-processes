function [pv, PartInfo, mRemaining, cRemaining] = ExchangePhases(...
    t, dt, iT, pv, mRemaining, cRemaining, PartInfo, SpeciesInfo, ModelParams, Const)
    %
    nPartGeneral = 2;
    iImmobile = 1;
    iMobileIni = 2;
    
    nT = numel(t);
    tOffset = t(iT);
    tAfter = t(iT:nT) - tOffset;
    tBounds = ModelParams.LogNorm.bounds;
    % Select only those time steps that are affected by current injection
    iCalcT = (tAfter <= tBounds(2));
    nCalcT = sum(iCalcT);
    iTend = iT + nCalcT - 1;
    tAfter = tAfter(iCalcT);

    mRemaining(1, iT + 1, SpeciesInfo.iInert) = mRemaining(1, iT, SpeciesInfo.iInert);
    mRemaining(2, iT + 1, :) = mRemaining(2, iT, :);
    % Prepare inputs for exchange equation
    iActivePart = [iImmobile, iMobileIni, nPartGeneral + (iT:iTend)];
%     cPartOutR = cat(2, cRemaining(1, iT, SpeciesInfo.iFlush), ...
%         PartInfo.GetConcentration(1, iActivePart, SpeciesInfo.iFlush));
%     pvPartOutR = cat(2, pv(1), PartInfo.GetVolume(1, iActivePart));
    cPartOutR = PartInfo.GetConcentration(1, iActivePart, SpeciesInfo.iFlush);
    pvPartOutR = PartInfo.GetVolume(1, iActivePart);
    % Don't calculate for particles vith (almost) zero volume
    iCalcLog = (pvPartOutR > Const.VOLUME_EPSILON);
    % ... but make sure immobile phase is calculated
    iCalcLog(1) = iImmobile;
    iCalc = find(iCalcLog);
    % Solve exchange equation in order to obtain concentrations in both phases at the end of
    % current time step. Different volumes will have different effect on exchange here
    if (isequal(iCalc, 1))
        cPartR = cPartOutR(:, iCalcLog, :);
    else
        cPartR = MultiConcentrationExchange([tAfter(1), tAfter(1) + dt], ...
            cPartOutR(:, iCalcLog, :), ModelParams.kExch, pvPartOutR(:, iCalcLog), Const);
    end
    % Update computed concentrations for particles
    PartInfo = PartInfo.SetConcentration(cPartR, ...
        1, iActivePart(iCalcLog), SpeciesInfo.iFlush);
    % Update total masses of solutes in phases
    mRemaining(1, iT + 1, SpeciesInfo.iFlush) = ...
        PartInfo.GetMass(1, iActivePart(1), SpeciesInfo.iFlush);
    mRemaining(2, iT + 1, SpeciesInfo.iFlush) = ...
        sum(PartInfo.GetMass(1, iActivePart(2:end), SpeciesInfo.iFlush), 2);

    return
    
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
        
        iIm = 1;
        iM = 2:nEl;
        Cim0 = cIni(1, iIm, :);
        vim = v(iIm);
        vm = sum(v(iM));
        Cm0 = nan(1, 1, nSolutes);
        for iSol = 1:nSolutes
            Cm0(1, 1, iSol) = sum(cIni(1, iM, iSol) .* v(iM)) / vm;
        end

        % Immobile phase
        t = tRange(2) - tRange(1);
        % Change in concentration of immobile phase
        Cim = (Cm0 * vm + Cim0 * vim) / (vm + vim) + vm * (-Cm0 + Cim0) * exp(-kExch * t) / (vm + vim);
        % Change in total concentration of all mobile elements
        Cm = (Cm0 * vm + Cim0 * vim) / (vm + vim) - vim * (-Cm0 + Cim0) * exp(-kExch * t) / (vm + vim);
        % Combine in the resulting array
        cPart = nan(1, nEl, nSolutes);
        cPart(1, iIm, :) = Cim;

        % Change of gradient between mobile-immobile concentrations
        % k = (Cm - Cim) ./ (Cm0 - Cim0);
        k = exp(-kExch * t) * ones(1, nSolutes);
        % Adjust individual mobile elements' concentrations to meet newly calculated mean
        for iSol = 1:nSolutes
            cPart(1, iM, iSol) = k(iSol) * (cIni(1, iM, iSol) - Cm0(iSol)) + Cm(iSol);
        end
    end
end