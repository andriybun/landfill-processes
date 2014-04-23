function MainMultiConcentrationExchangeBiochem
    % Testing method to calculate exchange of solute between different entities containing it. The
    % driving force of exchange is concentration gradient.
    % Dimensions of arrays: nElements x nTimeSteps x nSolutes
    %   nElements = (main phase + all particles it's exchanging with)

    clc;
    
    addpath('../../Common/');

    Const.EPSILON = 1e-9;

    dt = 1;
    tEnd = 50;
    tRange = 0:dt:tEnd;
    nT = numel(tRange);
    
    nEl = nT;
    pvAllEl = 1.5;
    
    nSpecies = 25;
    
    rngSeed = 1;
    rand(rngSeed);
    
    % Generate random rain data
%     isRainfall = (rand(1, tEnd) < 0.12);
%     rainfallIntensity = sum(reshape(isRainfall, dt, []), 1);
%     rainfallIntensity(rainfallIntensity < 1) = 0;
    rainfallIntensity = zeros(1, nT);
    rainfallIntensity(1:2) = 1;
%     rainfallIntensity(1:10:end) = 1;
%     rainfallIntensity(2:10:end) = 1;
    rainData = 5e-1 * rainfallIntensity;
    
    % Travel time pdf
    mu = 2;
    sigma = 0.7;
    tT = TtPdf(tRange, mu, sigma);

    % Initialize object to keep information about volumes and concentrations
    % dimensions of arrays (entry time x leave time)
    partVolume = zeros(nT, nT);

    % PartInfo = ConcentrationCl(zeros(1, nT), zeros(1, nT, nSpecies));
    % PartInfo = PartInfo.AddVolume(qOutAfter, :, iT:iTend);
    % PartInfo.GetVolume(1, iT:iTend);
    % PartInfo = PartInfo.SetConcentration(cPartR(end, 2:(nCalcT+1), :), ...
    %          1, iT:iTend, iFlushSpecies);
    % PartInfo.GetMass(1, iT, iFlushSpecies);
    
    % Numerical convolution
    qOut = zeros(size(rainData));
    for iT = 1:numel(tRange)
        t = tRange(iT);
        tAfter = tRange(iT:end) - t;
        if (rainData(iT) ~= 0)
            partVolume(iT, iT:end) = rainData(iT) * TtPdf(tAfter, mu, sigma);
            qOut(iT:end) = qOut(iT:end) + partVolume(iT, iT:end);
        end
    end

    % plotyy(tRange, rainData, tRange, qOut);
    
    cIni = zeros(1, nEl+1, nSpecies);
    cIni(1, 1, :) = 1;
    cIni(1, 2:(nEl+1), :) = repmat(rand(1, nEl), [1, 1, nSpecies]);
    
%     disp(squeeze(cIni));
    
    pv = zeros(1, nEl+1, 1);
    pv(1, 1) = 1;
    pv(1, 2:(nEl+1)) = partVolume(1, :);
%     rand(rngSeed);
%     pv(1, 2:(nEl+1)) = rand(1, nEl);
%     pv(1, 2:(nEl+1)) = pv(1, 2:(nEl+1)) * pvAllEl / sum(pv(1, 2:(nEl+1)));

    kExch = 1e-2;
    
    mIni = squeeze(sum(cIni .* repmat(pv, [1, 1, nSpecies]), 2));
    
    tic
    cPart = MultiConcentrationExchangeOdeBiochem(tRange, cIni, kExch, pv, partVolume, Const);
    toc
   
    mEnd = squeeze(sum(cPart(nT, :, :) .* repmat(pv, [1, 1, nSpecies]), 2));
    
%     iErr = reshape(find(abs(mIni - mEnd) >= Const.EPSILON), 1, []);
%     if (~isempty(iErr))
%         for iSolute = iErr
%             warning('\n\t%d-th solute: wrong mass balance. Error = %f', ...
%                 iSolute, abs(mIni(iSolute) - mEnd(iSolute)));
%         end
%     end
    
    vTot = zeros(1, nT);
    for iT = 1:nT
        vTot(iT) = sum(sum(partVolume(1:iT, iT+1:end)));
    end

    figH = figure(1);
    iSolute = 3;
    set(figH, 'Position', [1050, 200, 600, 450]);
    plot(tRange, squeeze(squeeze(cPart(sub2ind(size(cPart), 1:nT, 1 + (1:nT), repmat(iSolute, [1, nT]))))));
    hold on
    plot(tRange, squeeze(cPart(:, 1, iSolute)), 'r');
    hold off
    
    figure();
    plot(tRange, vTot);
%     plot(tRange, squeeze(squeeze(cPart(:, 2:end, iSolute))));
%     hold on
%     plot(tRange, squeeze(squeeze(cPart(:, 1, iSolute))), 'b', 'LineWidth', 2);
%     hold off
   
    return
    
    function ttX = TtPdf(tX, mu, sigma)
        ttCdf = logncdf(tX, mu, sigma);
        ttX = cat(2, 0, diff(ttCdf) ./ diff(tX));
    end
        
end