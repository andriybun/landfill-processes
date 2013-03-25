function DisplayPlotSinglePhase(timeOutVec, ...
                                cOutArr, ...
                                cAnalyticalArr, ...
                                ModelDim, ...
                                SoilPar, ...
                                inFlow, ...
                                selectedNodes, ...
                                doDisplayAnalyticalSolution)
    % Plot times series as a function of time for selected depths.
    % Compare analytical and numerical solutions.
    
    % Create figure with specified width (in order to fit legend)
    figH = figure(1);
    figPos = get(figH, 'Position');
    newWidth = 800;
    figPos(1) = figPos(1) + ceil((figPos(3) - newWidth) / 2);
    figPos(3) = newWidth;
    set(figH, 'Position', figPos);

    % Create custom color order matrix
    nCol = numel(selectedNodes);
    colorMatrix = zeros(nCol, 3);
    black2red = floor(nCol / 2);
    red2yellow = nCol - black2red;
    colorMatrix(1:black2red, 1) = linspace(0.2, 0.8, black2red);
    colorMatrix(black2red+1:nCol, 1) = 1;
    colorMatrix(black2red+1:nCol, 2) = linspace(0.1, 0.9, red2yellow);
    set(gca, 'ColorOrder', colorMatrix, 'nextplot', 'replacechildren');
    
    % Plot analytical solution
    if doDisplayAnalyticalSolution
        plot(timeOutVec, cAnalyticalArr(selectedNodes, :), 'LineWidth', 2);
    end

    % Plot numerical solution
    hold on
    soluteIdx = 1;
    plot(timeOutVec, squeeze(cOutArr(selectedNodes, soluteIdx, :)), 'o-', 'LineWidth', 1);
    hold off

    % Generate legend
    leg = cell(1, nCol);
    idx = 1;
    for pos = selectedNodes
        leg{idx} =          sprintf('%5.3f', ModelDim.zn(pos));
        idx = idx + 1;
    end
    legH = legend(leg, 'Location', 'EastOutside');
    legTitle = get(legH, 'Title');
    set(legTitle, 'string', 'Depth');
    
    % Adjust axes
    axBounds = axis;
    axis([axBounds(1:2), 0, 1]);
    xlabel('time [days]');
    ylabel('Concentration [m^3/m^3]');
    title(sprintf('Concentrations at depths (v = %5.3e, d = %5.3e)', inFlow, SoilPar.d(1)), 'FontWeight','bold');

end