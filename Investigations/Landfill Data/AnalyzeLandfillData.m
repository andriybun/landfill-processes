function AnalyzeLandfillData
    close all

    % Read data
    LeachData = load('data/WieringermeerProcessed');
    LevelData = load('data/WieringermeerLevelProcessed');
    RainData = load('data/RaindataProcessed');

    % Get initial time value. For convenience we count days from zero. Thus we need to know
    % initial time
    tBaseLeach = LeachData.data(1, 1);
    
    % Calculate daily values
    LeachDaily = DailyValues(LeachData, 'last');
    LevelDaily = DailyValues(LevelData, 'last');

    tic
    CombinedDaily = JoinArrays(LeachDaily, 1, LevelDaily, 1);
    CombinedDaily = JoinArrays(CombinedDaily, 1, RainData, 1);
    toc

    % Shift time and normalize some other data so that cumulative datasets start with zero
    for i = 1:6
        CombinedDaily.data(:, i) = CombinedDaily.data(:, i) - CombinedDaily.data(1, i);
    end

    CombinedDaily.data(:, 2) = [0; diff(CombinedDaily.data(:, 2))];
    
    plot(CombinedDaily.data(:, 1), 13 * cumsum(CombinedDaily.data(:, 2)), 'g');
    hold on
    plot(CombinedDaily.data(:, 1), cumsum(CombinedDaily.data(:, 9)), 'b');
    plot(CombinedDaily.data(:, 1), cumsum(CombinedDaily.data(:, 11)), 'r');
    legend('on site rain', 'meteo rain', 'meteo rain - evap', 'location', 'SouthEast');
    hold off
    
    iCell = 1;
    iRain = 11;
    
    close all
    figure(1);
    plotyy(CombinedDaily.data(:, 1), cumsum(CombinedDaily.data(:, iRain)), ...
        CombinedDaily.data(:, 1), CombinedDaily.data(:, 6 + iCell));
    legend('rain', 'level');
    figure(2);
    plotyy(CombinedDaily.data(:, 1), [cumsum(CombinedDaily.data(:, iRain)), CombinedDaily.data(:, 4 + iCell) * 0.1], ...
        CombinedDaily.data(:, 1), CombinedDaily.data(:, 6 + iCell));
    legend('rain', 'leachate pump', 'location', 'SouthEast');
    % Interval with constant level
    if iCell == 1
        iSel = 386:514;
    elseif iCell == 2
        iSel = 379:514;
    end
    figure(3);
    plotyy(CombinedDaily.data(iSel, 1), CombinedDaily.data(iSel, iRain), ...
        CombinedDaily.data(iSel, 1), [0; diff(CombinedDaily.data(iSel, 4 + iCell))]);
    title('Rainfall vs. leachate volume (constant level)');
    legend('rain', 'leachate');
    figure(4);
    plotyy(CombinedDaily.data(iSel, 1), cumsum(CombinedDaily.data(iSel, iRain)), ...
        CombinedDaily.data(iSel, 1), CombinedDaily.data(iSel, 4 + iCell));
    title('Rainfall vs. leachate volume (cumulative, constant level)');
    legend('rain', 'leachate');
    % Pump switched off
    if iCell == 1
        iSel = 165:325;
    elseif iCell == 2
        return
    end
    figure(5)
    plotyy(CombinedDaily.data(iSel, 1), CombinedDaily.data(iSel, iRain), ...
        CombinedDaily.data(iSel, 1), CombinedDaily.data(iSel, 6 + iCell));
    title('Rainfall vs. leachate level (pump off)');
    legend('rain', 'level');
    
    return
    
    % Some checking commands
    % Check percentage of NaN's (has to be as close to zero as possible): 
    %       sum(isnan(CombinedDaily.data(:, iRain))) / size(CombinedDaily.data(:, iRain), 1)
    % Display headers: 
    %       cat(2, num2cell(1:numel(CombinedDaily.headers))', CombinedDaily.headers)
    tDaily = CombinedDaily.data(:, 1);
    
    % Get required columns
    iDebit = [8];
    iLevel = iDebit + 4;
    qPumped = sum(diff(cat(1, CombinedDaily.data(:, iDebit), CombinedDaily.data(end, iDebit))), 2);
    lWell = smooth(CombinedDaily.data(:, iLevel), 1);
    dlWell = sum(diff(cat(1, lWell, lWell(end))), 2);
    
    % Translate level and pumped volume data into pumped volume
    % Coefficient to translate change in level of leachate to actual volume. This takes into account
    % total pore volume occupied by leachate inside the landfill based on its area and porosity.
    kVol = 0.04 * 500;
    qTot = kVol * dlWell + qPumped;
    
    
    % Plot volume fluxes vs. level data
    figure(1);
    plotyy(tDaily, [qPumped, qTot], tDaily, lWell);
    legend('Pumped', 'Total', 'Level');

%     % Plot net infiltration vs. leachate discharge
%     CombinedDailyAll = JoinArrays(CombinedDaily, 1, RainData, 1);
%     netInfiltration = CombinedDailyAll.data(:, 23);
%     figure(2); 
%     plotyy(tDaily, netInfiltration, tDaily, qTot);

%     % Plot rainfall and evapotranspiration data
%     figure(3);
%     plot(RainData.data(:, 1), RainData.data(:, 2:end));
%     legend(RainData.headers{2:end});
% 
%     % Plot Cumulative 
%     figure(4);
%     plot(RainData.data(:, 1), cumsum(RainData.data(:, 2:end), 1));
%     legend(RainData.headers{2:end});
    
%     % Deconvolution travel time PDF
%     sel = 1:100;
%     tPdf = ifft(fft(qTot(sel)) ./ fft(netInfiltration(sel)));
%     figure(5);
%     plot(tDaily(sel), tPdf);
    

end