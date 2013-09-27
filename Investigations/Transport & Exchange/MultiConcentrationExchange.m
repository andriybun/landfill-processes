function cPart = MultiConcentrationExchange(tauX, cIniX, kExchX, pvX, Const)
    % Function to do calculate exchange of solutes between main container and different 
    % particles residing in it. The first element along elements' axis is the main container.
    % All other particles exchange with it, tending towards concentration equilibrium.
    % Input parameters:
    %   tauX        - time span for calculations
    %   cIniX       - initial concentrations of solutes in main container (cIniX(1, 1, :)) and
    %                 other particles (cIniX(2:end, 1, :)). Dimensions of array are:
    %                 nTimeSteps x (nElements + 1) x nSolutes
    %   kExchX      - exchange coefficient between phases
    %   pvX         - volumes of main container (pv(1)) and particles (pv(2:end))

    % Get dimensions of a problem
    [~, nElX, nSolutes] = size(cIniX);
    nElX = nElX - 1;
    nTX = numel(tauX);

    % Initialize output array
    cPart = zeros(nTX, nElX+1, nSolutes);
    % First element along time axis is initial concentration
    cPart(1, :, :) = cIniX;

    % Find elements with zero volume.
    isPvZero = RealEq(pvX, 0, Const.EPSILON);
    isPvZero(1) = 0;
    isPvNotZero = ~isPvZero;
    isPvNotZero(1) = 0;

    if (sum(isPvNotZero) == 0)
        % If all elements have zero volume, they don't exist and we don't need to do 
        % calculations. Just repeat initial concentrations for all time steps
        cPart(2:nTX, 1, :) = repmat(cIniX(1, 1, :), [nTX-1, 1, 1]);
    else
        % Calculate total volume of all elements. This will be needed for computing average
        % concentration
        pvAll = sum(pvX(2:nElX+1));
        % And perform calculations for all time steps
        for iTX = 2:nTX
            tX = tauX(iTX);
            cPart(iTX, 1, :) = cPart(iTX-1, 1, :);
            mElTot = zeros(1, 1, nSolutes);
            dmEl = zeros(1, nElX, nSolutes);

            % First we run through all the elements and calculate total masses of solutes in
            % particles and potential mass exchange between them and main container given
            % difference in concentration
            for iEl = 1:nElX
                if isPvNotZero(iEl+1)
                    mElTot = mElTot + cPart(iTX-1, iEl+1, :) * pvX(iEl+1);
                    dmEl(1, iEl, :) = pvX(1) .* pvX(iEl+1) ./ (pvX(1) + pvX(iEl+1)) .* ...
                        (cPart(iTX-1, 1, :) - cPart(iTX-1, iEl+1, :)) .* ...
                        kExchX * (tX - tauX(iTX-1));
                end
            end
            % Then calculate mass exchange assuming average concentration of all the elements
            % with the main container
            dmAll = pvX(1) * pvAll / (pvX(1) + pvAll) * ...
                (cPart(iTX-1, 1, :) - mElTot / pvAll) .* kExchX .* (tX - tauX(iTX-1));
            dmAllSum = sum(dmEl, 2);
            for iSoluteX = 1:nSolutes
                % For each particle and solute mass exchange is adjusted to match averaged mass 
                % exchange, but keep relative mass change of particles constant 
                if ~RealEq(dmAll(1, 1, iSoluteX), 0, Const.EPSILON)
                    if RealEq(dmAllSum(1, 1, iSoluteX), 0, Const.EPSILON)
                        dmAdj = zeros(1, nElX, 1);
                    else
                        dmAdj = dmEl(1, :, iSoluteX) .* dmAll(1, 1, iSoluteX) ./ ...
                            dmAllSum(1, 1, iSoluteX);
                    end
                else
                    dmAdj = dmEl(1, :, iSoluteX);
                    dmAll(1, 1, iSoluteX) = dmAllSum(1, 1, iSoluteX);
                end
                % Update resulting concentrations based on computed before data
                cPart(iTX, isPvNotZero, iSoluteX) = cPart(iTX-1, isPvNotZero, iSoluteX) + ...
                    dmAdj(isPvNotZero(2:nElX+1)) ./ pvX(isPvNotZero);
                cPart(iTX, isPvZero, iSoluteX) = 0;
            end
            cPart(iTX, 1, :) = cPart(iTX, 1, :) - dmAll ./ pvX(1);
        end
    end
end