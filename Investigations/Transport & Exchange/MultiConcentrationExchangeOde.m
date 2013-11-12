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
    nTX = numel(tRange);

    cIni = permute(cIni, [2, 3, 1]);
    cIni = reshape(cIni, [], 1);
    [t, cResRaw] = ode45(@(tX, cX) dC(tX, cX, v', kExch, [1, nEl, nSolutes]), tRange, cIni);
    cPart = reshape(cResRaw, [nTX, nEl, nSolutes]);
    
    return
    
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
            % The main relationships
            dCdt(iSoluteOffset + 1) = sum(kExch * vX(iPartMob) ./ sumVx .* gradC - ...
                R(tX, cX(iSoluteOffset + 1)));
            dCdt(iPartMobSol) = -kExch * vX(1) ./ sumVx .* gradC;
        end
    end

    % Stub. Reactions function
    function rX = R(tX, cX)
        kX = 1e-3;
        rX = kX * (1 - cos(tX)) * cX;
    end
end