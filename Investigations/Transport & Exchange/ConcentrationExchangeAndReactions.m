function ConcentrationExchangeAndReactions
    % This is a script to test integration of exchange of concentrations with source/sink term
    % 

    % Constants
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
    
    % Other
    kExch = 1e-1;
    
    % Resulting arrays
    tic
    cRes = SolveOde(cIni, NUM_PHASES, nT, nSolutes, kExch);
    toc
    
    for iSolute = 1:nSolutes
        subplot(2, 2, iSolute);
        plot(t, squeeze(cRes(:, :, iSolute)));
    end
    
    return
    
    % Stub. Reactions function
    function rX = R(tX, cX)
        kX = 1e-1;
        rX = kX * (1 - cos(tX)) * cX;
    end

    % Wrapper function
    function cRes = SolveOde(cIni, NUM_PHASES, nT, nSolutes, kExch)
        cIniX = permute(cIni, [1, 3, 2]);
        cIniX = reshape(cIniX, [], 1);
        [t, cResRaw] = ode45(@(tX, cX) dC(tX, cX, kExch), tRange, cIniX);
        cRes = permute(reshape(cResRaw, [nT, NUM_PHASES, nSolutes]), [2, 1, 3]);
    end

    % System of ODE's
    function dCdt = dC(tX, cX, kExch)
        nEl = numel(cX);
        dCdt = zeros(nEl, 1);
        mobEl = 2:2:nEl;
        immobEl = 1:2:nEl;
        gradC = cX(mobEl) - cX(immobEl);
        dCdt(immobEl) = kExch * gradC - R(tX, cX(immobEl));
        dCdt(mobEl) = -kExch * gradC;
    end
end