function cPart = MultiConcentrationExchangeOdeBiochem(tRange, cIni, kExch, v, partVolume, Const)
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
    [~, nEl, ~] = size(cIni);
    nT = numel(tRange);
    
    %% Biochem ini
    % Initializing chemical module
    CHEM_MODEL_DIR = '../';
    addpath(genpath(CHEM_MODEL_DIR));
    [mIni, Comp, Pm, S, Rp] = initialize_ODE([CHEM_MODEL_DIR 'Pmatrix/Pmatrix.csv']);
%     ORI = initialize_ORI(Comp, 0);
    
    mInertIni = [1];
    nInertSolutes = numel(mInertIni);
    nReactiveSolutes = numel(mIni);
    mIni = cat(2, mIni, mInertIni);
    nSolutes = numel(mIni);
    iReactiveSolutes = 1:nReactiveSolutes;
    iInertSolutes = nReactiveSolutes:nSolutes;
    iFlushSolutes = [2:5, 8:9, 22];
    nFlushSolutes = numel(iFlushSolutes);
    
    % Initial concentrations of solutes
    vBiochem = 0.325;
    vAdj = v / vBiochem;
    
    mIni = cat(2, permute(mIni, [1, 3, 2]) * vAdj(1), zeros(1, nEl-1, nSolutes));
    cIni = mIni ./ repmat(v, [1, 1, nSolutes]);
    cIni(:, RealEq(v, 0, Const.EPSILON), :) = 0;
    %% Biochem ini end
    
    global Call V2 tt Rall
    Call = [];  V2 = []; tt = []; Rall = [];
    options = odeset('OutputFcn', ...
        'Store_Orchestra_Results', ...
        'AbsTol', 1e-6, ...
        'RelTol', 1e-6);
%         'Refine', 1, ...

    % Integrate volumes of particles before ODE solver
    vSum = zeros(nT, nT);
    
    for iT = 1:nT
        if iT == 1
            vSum(iT, iT:nT) = partVolume(iT, iT:nT);
        else
            vSum(iT, iT:nT) = vSum(iT-1, iT:nT) + partVolume(iT, iT:nT);
        end
    end
    partVolume = cat(1, partVolume, zeros(1, nT));
    
    cIni = permute(cIni, [2, 3, 1]);
    cIni = reshape(cIni, [], 1);
    [t, cResRaw] = ode15s(@(tX, cX) dC(tX, cX, kExch, [1, nEl, nSolutes]), [tRange(1), tRange(end)], cIni, options);
    
    cPart = reshape(cResRaw, [nT, nEl, nSolutes]);
    
    return
    
    % System of ODE's
    function dCdt = dC(tX, cX, kExch, dimVec)
        nEl = dimVec(2);
        nSolutes = dimVec(3);
        iT = find(tX >= tRange, true, 'last');
        
        cX = cX';
        vX = cat(2, v(1), vSum(iT, :));
        
        % Resulting array
        dCdt = zeros(nEl * nSolutes, 1);
        % Indices of mobile particles
        iPartMob = 2:nEl;
        % Loop over solutes (this is slightly faster than matrix operations
        for iSoluteR = 1:nFlushSolutes
            iSolute = iFlushSolutes(iSoluteR);
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
            % Immobile concentration
            dCdt(iSoluteOffset + 1) = sum(kExch * vX(iPartMob) ./ sumVx .* gradC);
            % Exclude particles with zero volume
            doCalc = ~RealEq(vSum(iT, :), 0, Const.EPSILON);
            % Mobile concentrations
            Q = cX(iPartMobSol(doCalc)) .* (-partVolume(iT+1, doCalc) ./ ...
                (vSum(iT, doCalc) + partVolume(iT+1, doCalc)));
            dCdt(iPartMobSol(doCalc)) = -kExch * vX(1) ./ sumVx(doCalc) .* gradC(doCalc) + Q;
        end
        iPartImmob = ((1:nSolutes) - 1) * nEl + 1;
        iPartImmob = iPartImmob(iReactiveSolutes');
        dCdt(iPartImmob) = dCdt(iPartImmob) + R(tX, cX(iPartImmob), vX)';
    end

    % Stub. Reactions function
    function rX = R(tX, cX, vX)
%         kX = 1e-3;
%         rX = kX * (1 - cos(tX)) * cX;
        rX = cX * 0;
%         mX = cX * vBiochem;
%         dM = bioreactor(tX, mX, Comp, Pm, S, Rp, ORI);
%         rX = dM ./ vX(1);
    end
end