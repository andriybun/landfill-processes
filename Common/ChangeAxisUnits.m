function ChangeAxisUnits(hFig, axisName, unitScale)

    axesObjs = get(hFig, 'Children');
    dataObjs = get(axesObjs, 'Children');
    for idx = 2:numel(dataObjs)
        dataObj = dataObjs{idx};
        axisLimStr = sprintf('%slim', axisName);
        oldLims = get(axesObjs(idx), axisLimStr) * unitScale;
        set(axesObjs(idx), axisLimStr, [-Inf, Inf]);
%         newLims = get(axesObjs(idx), axisLimStr);
%         set(axesObjs(idx), axisLimStr, [max(oldLims(1), newLims(1)), min(oldLims(2), newLims(2))]);
        if any(strcmp(get(dataObj, 'Type'), 'line'))
            axisData = get(dataObj, sprintf('%sData', axisName));
            if iscell(axisData)
                for dataIdx = 1:numel(axisData)
                    axisDataOut = axisData{dataIdx} * unitScale;
                    set(dataObj(dataIdx), sprintf('%sData', axisName), axisDataOut);
                end
            else
                axisDataOut = axisData * unitScale;
                set(dataObj, sprintf('%sData', axisName), axisDataOut);
            end
        end
    end
    
end
