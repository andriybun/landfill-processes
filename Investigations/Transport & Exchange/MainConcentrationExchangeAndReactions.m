function MainConcentrationExchangeAndReactions
    % This is a script to test integration of exchange of concentrations with source/sink term
    % 

    % Constants
    EPSILON = 1e-10;
    NUM_PHASES = 2;
    
    % Time parameters
    tEnd = 100;
    dt = 0.1;
    tRange = 0:dt:tEnd;
    nT = numel(tRange);
    
    % Initial concentrations (we distinguish 2 phases: 1 - immobile, 2 - mobile)
    nSolutes = 3;
    cIni = [0, 0.9, 1;
            1, 0.1, 1];
    cIni = permute(cIni, [1, 3, 2]);
    
    % Volumes of phases
    v = [9; 1];
    
    % Other
    kExch = 1e-1;
    
    % Resulting arrays
    tic
    cRes = SolveOde(cIni, v, NUM_PHASES, nT, nSolutes, kExch);
    toc
    
    % Check mass balance
    mIni = sum(Mass(v, cIni), 1);
    mEnd = sum(Mass(v, cRes(:, end, :)), 1);
    fprintf('Maximum mass balance error is %e\n', max(abs(reshape(mEnd - mIni, [], 1))));
    
    % Plot
    for iSolute = 1:nSolutes
        subplot(2, 2, iSolute);
        plot(t, squeeze(cRes(:, :, iSolute)));
    end
    
    return
    
    % Wrapper function
    function cRes = SolveOde(cIni, v, NUM_PHASES, nT, nSolutes, kExch)
        cIniX = permute(cIni, [1, 3, 2]);
        cIniX = reshape(cIniX, [], 1);
        [t, cResRaw] = ode45(@(tX, cX) dC(tX, cX, v, kExch), tRange, cIniX);
        cRes = permute(reshape(cResRaw, [nT, NUM_PHASES, nSolutes]), [2, 1, 3]);
    end

    % System of ODE's
    function dCdt = dC(tX, cX, vX, kExch)
        nEl = numel(cX);
        dCdt = zeros(nEl, 1);
        mobEl = 2:2:nEl;
        immobEl = 1:2:nEl;
        gradC = cX(mobEl) - cX(immobEl);
        sumVx = vX(1) + vX(2);
        dCdt(immobEl) = kExch * vX(2) / sumVx * gradC - R(tX, cX(immobEl));
        dCdt(mobEl) = -kExch * vX(1) / sumVx * gradC;
    end

    % Stub. Reactions function
    function rX = R(tX, cX)
        kX = 1e-1;
        rX = kX * (1 - cos(tX)) * cX;
    end

    % Compute masses of solutes
    function mX = Mass(vX, cX)
        [nPhasesX, ~, nSolutesX] = size(cX);
        mX = nan(nPhasesX, 1, nSolutesX);
        for iSoluteX = 1:nSolutesX
            mX(:, 1, iSoluteX) = vX .* cX(:, 1, iSoluteX);
        end
    end

end