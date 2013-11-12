function MainMultiConcentrationExchange
    % Testing method to calculate exchange of solute between different entities containing it. The
    % driving force of exchange is concentration gradient.
    % Dimensions of arrays: nElements x nTimeSteps x nSolutes
    %   nElements = (main phase + all particles it's exchanging with)

    clc;
    
    addpath('../../Common/');

    Const.EPSILON = 1e-9;

    dt = 0.01;
    tEnd = 10;
    tRange = 0:dt:tEnd;
    nT = numel(tRange);
    
    nEl = 500;
    pvAllEl = 4;
    
    nSolutes = 2;
    
%     rngSeed = 1;
%     rand(rngSeed);
    
    cIni = zeros(1, nEl+1, nSolutes);
    cIni(1, 1, :) = 1;
    cIni(1, 2:(nEl+1), :) = repmat(rand(1, nEl), [1, 1, nSolutes]);
    
%     disp(squeeze(cIni));
    
    pv = zeros(1, nEl+1, 1);
    pv(1, 1) = 1;
%     rand(rngSeed);
    pv(1, 2:(nEl+1)) = rand(1, nEl);
    pv(1, 2:(nEl+1)) = pv(1, 2:(nEl+1)) * pvAllEl / sum(pv(1, 2:(nEl+1)));

    kExch = 2;
    
    mIni = squeeze(sum(cIni .* repmat(pv, [1, 1, nSolutes]), 2));
    
    tic
    cPart = MultiConcentrationExchangeOde(tRange, cIni, kExch, pv, Const);
    toc
    
    mEnd = squeeze(sum(cPart(nT, :, :) .* repmat(pv, [1, 1, nSolutes]), 2));
    
    iErr = reshape(find(abs(mIni - mEnd) >= Const.EPSILON), 1, []);
    if (~isempty(iErr))
        for iSolute = iErr
            warning('\n\t%d-th solute: wrong mass balance. Error = %f', ...
                iSolute, abs(mIni(iSolute) - mEnd(iSolute)));
        end
    end
    
    figH = figure(1);
    set(figH, 'Position', [1050, 200, 600, 450]);
    plot(tRange, squeeze(squeeze(cPart(:, 2:end, 1))));
    hold on
    plot(tRange, squeeze(squeeze(cPart(:, 1, 1))), 'b', 'LineWidth', 2);
    hold off
    return
    
end