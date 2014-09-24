function cPart = MultiConcentrationExchangeOde(tRange, cIni, kExch, v, Const)
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
    % Total change of mass of compounds in all mobile elements
    dm = (Cm - Cm0) .* vm;
    % Distribute this change over all elements proportionally to 1/vm
    dmDistr = nan(1, nEl, nSolutes);
    % Combine in the resulting array
    cPart = nan(1, nEl, nSolutes);
    cPart(1, 1, :) = Cim;
    
    % Change of gradient between mobile-immobile concentrations
    k = (Cm - Cim) ./ (Cm0 - Cim0);
    % Adjust individual mobile elements' concentrations to meet newly calculated mean
    for iSol = 1:nSolutes
        cPart(1, iM, iSol) = k(iSol) * (cIni(1, iM, iSol) - Cm0(iSol)) + Cm(iSol);
    end
end