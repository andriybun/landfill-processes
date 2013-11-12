function Res = DailyValues(Data, STAT_TYPE)
    % Data      - structure with fields 'data' and 'headers' containing data
    % tIni      - initial time value that defines beginning of calculations
    % STAT_TYPE - type of statistics to be calculated for every day. May be one of the following:
    %             'LAST', 'SUM', 'AVG', 'MAX', 'MIN'
    
    % Remove NaN's
    dataNoNan = RemoveNan(Data.data);
    
    [nRowsOrig, nCols] = size(dataNoNan);
    
    % Shifted time vector
    t = dataNoNan(:, 1);
    
    % Coarse time discretization
    tCoarse = floor(t) + 1;
    % Unique time moments
    tDaily = unique(tCoarse);
    
    % Number of unique time moments
    nRowsDaily = numel(tDaily) - 1;
    
    % Initialize resulting structure
    Res = struct();
    Res.headers = Data.headers;
    Res.data = nan(nRowsDaily, nCols);
    
    switch upper(STAT_TYPE)
        case 'LAST'
            iDaily = 1;
            for i = 1:nRowsOrig-1
                if tCoarse(i) ~= tCoarse(i+1)
                    Res.data(iDaily, 1) = tCoarse(i);
                    Res.data(iDaily, 2:end) = dataNoNan(i, 2:end);
                    iDaily = iDaily + 1;
                end
            end
%         case 'SUM'
%             
%         otherwise
%             
    end
          
    return
    
    function resX = RemoveNan(inX)
        isNan = isnan(inX(:, 2));
        resX = inX(~isNan, :);
    end
end