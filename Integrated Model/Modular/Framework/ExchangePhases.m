function [pv, MobileC, LongTermC, ImmobileC, mRemaining, cRemaining] = ExchangePhases(...
    t, dt, iT, pv, mRemaining, cRemaining, MobileC, LongTermC, ImmobileC, ...
    SpeciesInfo, ModelParams, Const)

    nT = numel(t);
    tOffset = t(iT);
    tAfter = t(iT:nT) - tOffset;
    tBounds = ModelParams.LogNorm.bounds;
    % Select only those time steps that are affected by current injection
    iCalcT = (tAfter <= tBounds(2));
    nCalcT = sum(iCalcT);
    iTend = iT + nCalcT - 1;
    tAfter = tAfter(iCalcT);

    % Prepare inputs for exchange equation
%     cMob = MobileC.GetConcentration(1:iT, iT:iTend, SpeciesInfo.iFlush);
    cMob = MobileC.GetConcentration(1, iT:iTend, SpeciesInfo.iFlush);
    cLongTerm = LongTermC.GetConcentration(:, :, SpeciesInfo.iFlush);
    cImmob = ImmobileC.GetConcentration(:, :, SpeciesInfo.iFlush);
%     pvMob = MobileC.GetVolume(1:iT, iT:iTend);
    pvMob = MobileC.GetVolume(1, iT:iTend);
    pvLongTerm = LongTermC.GetVolume();
    pvImmob = ImmobileC.GetVolume();
    % Solve exchange equation in order to obtain concentrations in both phases at the end of
    % current time step. Different volumes will have different effect on exchange here
    [cMob, cLongTerm, cImmob] = MultiConcentrationExchange([tAfter(1), tAfter(1) + dt], ...
        cMob, cLongTerm, cImmob, ModelParams.kExch, pvMob, pvLongTerm, pvImmob);
    % Update computed concentrations for particles
%     MobileC = MobileC.SetConcentration(cMob, 1:iT, iT:iTend, SpeciesInfo.iFlush);
    MobileC = MobileC.SetConcentration(cMob, 1, iT:iTend, SpeciesInfo.iFlush);
    LongTermC = LongTermC.SetConcentration(cLongTerm, :, :, SpeciesInfo.iFlush);
    ImmobileC = ImmobileC.SetConcentration(cImmob, :, :, SpeciesInfo.iFlush);
    % Update total masses of solutes in phases
    mRemaining(1, iT + 1, SpeciesInfo.iFlush) = ...
        ImmobileC.GetMass(:, :, SpeciesInfo.iFlush);
%     mRemaining(2, iT + 1, SpeciesInfo.iFlush) = ...
%         LongTermC.GetMass(:, :, SpeciesInfo.iFlush) + ...
%         sum(sum(MobileC.GetMass(1:iT, iT:iTend, SpeciesInfo.iFlush), 2), 1);
    mRemaining(2, iT + 1, SpeciesInfo.iFlush) = ...
        LongTermC.GetMass(:, :, SpeciesInfo.iFlush) + ...
        sum(sum(MobileC.GetMass(1, iT:iTend, SpeciesInfo.iFlush), 2), 1);

    return
    
    function [cMob, cLongTerm, cImmob] = ...
            MultiConcentrationExchange(tRange, cMob, cLongTerm, cImmob, kExch, ...
            pvMob, pvLongTerm, pvImmob)
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
        [~, ~, nSolutes] = size(cMob);
        
        Cim0 = cImmob;
        vim = pvImmob;
        vm = sum(sum(pvMob)) + pvLongTerm;
        Cm0 = nan(1, 1, nSolutes);
        for iSol = 1:nSolutes
            Cm0(1, 1, iSol) = (sum(sum(cMob(:, :, iSol) .* pvMob, 2), 1) + ...
                cLongTerm(1, 1, iSol) * pvLongTerm) / vm;
        end

        % Immobile phase
        t = tRange(2) - tRange(1);
        % Change in concentration of immobile phase
        Cim = (Cm0 * vm + Cim0 * vim) / (vm + vim) + vm * (-Cm0 + Cim0) * exp(-kExch * t) / (vm + vim);
        % Change in total concentration of all mobile elements
        Cm = (Cm0 * vm + Cim0 * vim) / (vm + vim) - vim * (-Cm0 + Cim0) * exp(-kExch * t) / (vm + vim);
        % Results
        cImmob = Cim;

        % Change of gradient between mobile-immobile concentrations
        % k = (Cm - Cim) ./ (Cm0 - Cim0);
        k = exp(-kExch * t);
        % Adjust individual mobile elements' concentrations to meet newly calculated mean
        for iSol = 1:nSolutes
            cMob(:, :, iSol) = k * (cMob(:, :, iSol) - Cm0(iSol)) + Cm(iSol);
        end
        cLongTerm = k * (cLongTerm - Cm0) + Cm;
    end
end