function RunWithGui(appName, func, varargin)
% This function is used for running generic function using input parameters defined. It is also able
% to plot the result.
% Inputs:
%   appName - string. The GUI window will have the name as specified by this parameter.
%   func - handle of a function to run. The function may accept the inputs passed as 'input' 
%           structure and a progress bar [0, 100]. Its values can be set using 'pvalue' property.
%           The function must return a structure with results. Its fields then can be plotted using
%           GUI.
%   varargin - structures with parameters, that will be passed to the function. Note, composite 
%           types of fields (vectors, matrices, cells, sub-structures) may not work properly. In 
%           any case, it will be not possible to edit values of them.
% This function generates GUI window that allows change values for properties contained in input
% structure.
% After the function defined by handle finishes, a user can select its output data to be plot. For
% this the function should return a structure with vectors for each data entity.
    
    inputNames = cell(1, nargin-2);
    for iVar = 1:nargin-2
        inputNames{iVar} = inputname(iVar + 2);
    end
    inputCell = struct2cellX(varargin{:}, inputNames);
    CreateGui(inputCell, func);
    
    return
    
    function CreateGui(inputCell, func)
        % Window properties
        WINDOW_WIDTH = 400;
        WINDOW_HEIGHT = 630;
        % Progress bar properties
        PRBAR_HEIGHT = 20;
        PRBAR_MARGIN = 10;
        % Table properties
        MAX_TABLE_HEIGHT = 200;
        COLUMN_WIDTHS_REL = [90, 140, 70];
        ROW_HEIGHT = 19;
        BOUND_PIX = 2;
        % Button propeties
        BUTTON_WIDTH = 80;
        BUTTON_HEIGHT = 25;
        % Pop-up properties
        LIST_HEIGHT = 20;
        TEXT_HEIGHT = 13;
        TEXT_OFFSET = 5;
        
        % Create and then hide the GUI as it is being constructed.
        guiWindow = figure('Visible', 'off', 'Position', [200, 200, WINDOW_WIDTH, WINDOW_HEIGHT]);
        set(guiWindow, 'Name', appName, 'NumberTitle', 'off', 'Tag', 'guiWindow');
        set(guiWindow, 'Resize', 'off');
        windowColor = get(guiWindow, 'Color');
        
        % Progress bar
        prBar = progressbar(guiWindow, ...
            [PRBAR_MARGIN, PRBAR_MARGIN, WINDOW_WIDTH - 2 * PRBAR_MARGIN, PRBAR_HEIGHT], [0 100]);
        
        % Construct components
        % Button to run model
        btnVertPos = WINDOW_HEIGHT - PRBAR_MARGIN - BUTTON_HEIGHT;
        
        % Table
        tableWidth = WINDOW_WIDTH - 2 * PRBAR_MARGIN;
        columnWidths = round(tableWidth / sum(COLUMN_WIDTHS_REL) * COLUMN_WIDTHS_REL);
        columnWidths = columnWidths - [BOUND_PIX, sum(columnWidths) - tableWidth - 1, BOUND_PIX];
        if (MAX_TABLE_HEIGHT > ((size(inputCell, 1) + 1) * ROW_HEIGHT + BOUND_PIX * 2 - 1))
            tableHeight = (size(inputCell, 1) + 1) * ROW_HEIGHT + BOUND_PIX * 2 - 1;
        else
            tableHeight = MAX_TABLE_HEIGHT;
            columnWidths = columnWidths - [5, 5, 5];
        end
        tableVertPos = btnVertPos - PRBAR_MARGIN - tableHeight;
        hTable = uitable('Position', ...
            [PRBAR_MARGIN, tableVertPos, tableWidth, tableHeight]);
        set(hTable, 'ColumnName', {'Parameter', 'Value', 'Type'});
        set(hTable, 'ColumnEditable', [false, true, false]);
        set(hTable, 'ColumnWidth', num2cell(columnWidths));
        set(hTable, 'RowName', []);
        set(hTable, 'Data', inputCell(:, 1:3));
        
        % Pop-up lists
        popListTitlesVertPos = tableVertPos - PRBAR_MARGIN - TEXT_HEIGHT;
        popListVertPos = popListTitlesVertPos - TEXT_OFFSET - LIST_HEIGHT;
        popListWidth = (WINDOW_WIDTH - 3 * PRBAR_MARGIN) * 0.4;
        uicontrol('Style', 'text', ...
            'String', 'X-axis dataset',...
            'Backgroundcolor', windowColor, ...
            'Position', [PRBAR_MARGIN, popListTitlesVertPos, popListWidth, TEXT_HEIGHT]);
        hPopListX = uicontrol('Style', 'popupmenu', ...
            'Position', [PRBAR_MARGIN, popListVertPos, popListWidth, LIST_HEIGHT]);
        popListYhorOffset = 2 * PRBAR_MARGIN + popListWidth;
        uicontrol('Style', 'text', ...
            'String', 'Y-axis dataset',...
            'Backgroundcolor', windowColor, ...
            'Position', [popListYhorOffset, popListTitlesVertPos, popListWidth, TEXT_HEIGHT]);
        hPopListY = uicontrol('Style', 'popupmenu', ...
            'Position', [popListYhorOffset, popListVertPos, popListWidth, LIST_HEIGHT]);
        popList3rdDimHorOffset = 3 * PRBAR_MARGIN + 2 * popListWidth;
        popList3rdDimWidth = WINDOW_WIDTH - 4 * PRBAR_MARGIN - 2 * popListWidth;
        uicontrol('Style', 'text', ...
            'String', '3rd dim',...
            'Backgroundcolor', windowColor, ...
            'Position', [popList3rdDimHorOffset, popListTitlesVertPos, ...
            popList3rdDimWidth, TEXT_HEIGHT]);
        hPopList3rdDim = uicontrol('Style', 'popupmenu',...
            'Position', [popList3rdDimHorOffset, popListVertPos, popList3rdDimWidth, LIST_HEIGHT]);
        
        % Plot
        plotWidth = WINDOW_WIDTH - 60;
        plotHorPos = WINDOW_WIDTH - PRBAR_MARGIN - plotWidth;
        plotVertPos = 60 + PRBAR_HEIGHT;
        plotHeight = popListVertPos - plotVertPos - 2 * PRBAR_MARGIN;
        hPlotAxes = axes(...
            'Parent', guiWindow, ...
            'Units', 'pixels', ...
            'HandleVisibility', 'callback', ...
            'Position', [plotHorPos, plotVertPos, plotWidth, plotHeight]);
        
        % Add button to save data and figure to file
        btnSaveDataHorPos = 2 * PRBAR_MARGIN + BUTTON_WIDTH; 
        hBtnSaveData = uicontrol('Style', 'pushbutton', ...
            'String', 'Data to file', ...
            'Tag', 'hBtnSaveData', ...
            'Enable', 'off', ...
            'Position', [btnSaveDataHorPos , btnVertPos, BUTTON_WIDTH, BUTTON_HEIGHT]);
        btnSaveFigHorPos = 3 * PRBAR_MARGIN + 2 * BUTTON_WIDTH; 
        hBtnSaveFig = uicontrol('Style', 'pushbutton', ...
            'String', 'Plot to file', ...
            'Tag', 'hBtnSaveFig', ...
            'Enable', 'off', ...
            'Position', [btnSaveFigHorPos, btnVertPos, BUTTON_WIDTH, BUTTON_HEIGHT]);
        
        uiElements = struct();
        uiElements.guiWindow = guiWindow;
        uiElements.prBar = prBar;
        uiElements.hPopListX = hPopListX;
        uiElements.hPopListY = hPopListY;
        uiElements.hPopList3rdDim = hPopList3rdDim;
        uiElements.hPlotAxes = hPlotAxes;
        uiElements.hTable = hTable;
        uiElements.hBtnSaveData = hBtnSaveData;
        uiElements.hBtnSaveFig = hBtnSaveFig;
        
        ResetPopups(uiElements);
        
        set(hBtnSaveFig, 'Callback', {@SaveGraph, uiElements});
        set(hBtnSaveData, 'Callback', {@SaveData, uiElements});
        set(hPopListX, 'Callback', {@PlotGraph, uiElements});
        set(hPopListY, 'Callback', {@PlotGraph, uiElements});
        set(hPopList3rdDim, 'Callback', {@PlotGraph3d, uiElements});
        
        % Add button to run function
        hBtnStart = uicontrol('Style', 'pushbutton', ...
            'String', 'Start', ...
            'Tag', 'hBtnStart', ...
            'Position', [PRBAR_MARGIN, btnVertPos, BUTTON_WIDTH, BUTTON_HEIGHT], ...
            'Callback', {@FuncWrap, {func, inputCell, uiElements}});
        
        % Make the GUI visible.
        set(guiWindow, 'Visible', 'on');
    end
    
    function FuncWrap(hObject, eventdata, args)
        % Wrapper is used to run function with given set of arguments
        func = args{1};
        origInputs = args{2};
        uiElements = args{3};
        DeactivateGraph(uiElements);
        ResetPopups(uiElements);
        varArgPass = cell2structX(get(uiElements.hTable, 'data'), origInputs);
        output = func(varArgPass{:}, uiElements.prBar);
        set(uiElements.hPopListX, 'String', fieldnames(output));
        set(uiElements.hPopListX, 'Enable', 'on');
        set(uiElements.hPopListY, 'String', fieldnames(output));
        set(uiElements.hPopListY, 'Enable', 'on');
        set(hObject, 'UserData', output);
        cla(uiElements.hPlotAxes);
        set(uiElements.hBtnSaveData, 'Enable', 'on');
        guidata(hObject);
    end

    function PlotGraph(hObject, eventdata, uiElements)
        DeactivateGraph(uiElements);
        output = get(findobj('Tag', 'hBtnStart'), 'UserData');
        xStrings = get(uiElements.hPopListX, 'String');
        xAxisData = xStrings{get(uiElements.hPopListX, 'Value')};
        yStrings = get(uiElements.hPopListY, 'String');
        yAxisData = yStrings{get(uiElements.hPopListY, 'Value')};
        if ~isvector(output.(xAxisData))
            %% TODO: is error
        end
        switch (ndims(output.(yAxisData)))
            case 2
                plot(uiElements.hPlotAxes, output.(xAxisData), output.(yAxisData));
                ActivateGraph(uiElements);
            case 3
                set(uiElements.hPopList3rdDim, 'Enable', 'on');
                set(uiElements.hPopList3rdDim, 'String', ...
                    cat(2, {'all'}, num2cellStr(1:size(output.(yAxisData), 3))));
                set(uiElements.hPopList3rdDim, 'Value', 1);
        end
        xlabel(uiElements.hPlotAxes, xAxisData);
        ylabel(uiElements.hPlotAxes, yAxisData);
    end

    function PlotGraph3d(hObject, eventdata, uiElements)
        output = get(findobj('Tag', 'hBtnStart'), 'UserData');
        xStrings = get(uiElements.hPopListX, 'String');
        xAxisData = xStrings{get(uiElements.hPopListX, 'Value')};
        yStrings = get(uiElements.hPopListY, 'String');
        yAxisData = yStrings{get(uiElements.hPopListY, 'Value')};
        yDimStrings = get(uiElements.hPopList3rdDim, 'String');
        if strcmp(yDimStrings{get(uiElements.hPopList3rdDim, 'Value')}, 'all')
            [x, y] = meshgrid(output.(xAxisData), 1:size(output.(yAxisData), 3));
            mesh(uiElements.hPlotAxes, x, y, permute(output.(yAxisData), [3, 2, 1]));
        else
            yDimData = str2double(yDimStrings{get(uiElements.hPopList3rdDim, 'Value')});
            plot(uiElements.hPlotAxes, output.(xAxisData), ...
                squeeze(output.(yAxisData)(:, :, yDimData)));
        end
        ActivateGraph(uiElements);
    end

    function SaveData(hObject, eventdata, uiElements)
        [fileName, pathName] = uiputfile('*.mat', 'Select file to write data');
        if ~isequal(fileName, 0) 
            fileName = fullfile(pathName, fileName);
            data = get(findobj('Tag', 'hBtnStart'), 'UserData');
            save(fileName, '-struct', 'data');
        end
    end

    function SaveGraph(hObject, eventdata, uiElements)
        [fileName, pathName] = uiputfile('*.fig', 'Select file to write plot');
        if ~isequal(fileName, 0) 
            fileName = fullfile(pathName, fileName);
            SaveFigureToFile(uiElements.hPlotAxes, fileName);
        end
    end

    function ActivateGraph(uiElements)
        set(uiElements.hBtnSaveFig, 'Enable', 'on');
        drawnow;
    end

    function DeactivateGraph(uiElements)
        cla(uiElements.hPlotAxes);
        set(uiElements.hPopList3rdDim, 'Enable', 'off');
        set(uiElements.hBtnSaveFig, 'Enable', 'off');
        set(uiElements.hBtnSaveData, 'Enable', 'off');
        drawnow;
    end

    function ResetPopups(uiElements)
        set(uiElements.hPopListX, 'String', 'x-axis data', 'Enable', 'off', 'Value', 1);
        set(uiElements.hPopListY, 'String', 'y-axis data', 'Enable', 'off', 'Value', 1); 
        set(uiElements.hPopList3rdDim, 'String', '3rd dim', 'Enable', 'off', 'Value', 1);
        drawnow;
    end

    function outCell = struct2cellX(varargin)
        % Convert input structures into cell array of fields and their values
        inStruct = varargin(1:end-1);
        structName = varargin{end};
        % Initialize output struct
        outCell = cell(0, 4);
        % Process all structures
        for idx = 1:numel(inStruct)
            % Get fields of structure
            fieldNameVec = fieldnames(inStruct{idx});
            nFields = numel(fieldNameVec);
            % Append row with the name of struct:
            outCell = cat(1, outCell, {['# ' structName{idx}], '', 'input struct', ''});
            rowOffset = size(outCell, 1);
            % Append blank cells for data
            outCell = cat(1, outCell, cell(nFields, 4));
            % Process data
            for i = 1:nFields
                fieldName = fieldNameVec{i};
                fieldVal = inStruct{idx}.(fieldName);
                sz = size(fieldVal);
                if ~isnumeric(fieldVal) && ~ischar(fieldVal) && ~islogical(fieldVal)
                    fieldType = class(fieldVal);
                    valStr = sprintf('%s %d x %d', fieldType, sz(1), sz(2));
                    outCell(rowOffset + i, :) = ...
                        {fieldName, valStr, sprintf('#%s', fieldType), fieldVal};
                else
                    if isequal(sz, [1, 1]) || ischar(fieldVal)
                        outCell(rowOffset + i, 1:3) = {fieldName, fieldVal, class(fieldVal)};
                    else
                        valStr = sprintf('vector %d x %d', sz(1), sz(2));
                        classStr = sprintf('#%s', class(fieldVal));
                        outCell(rowOffset + i, :) = {fieldName, valStr, classStr, fieldVal};
                    end
                end
            end
        end
    end

    function outStruct = cell2structX(inCell, origCell)
        % Convert array of cells containing field names and their values into structures to be
        % passed to function
        outStruct = cell(1, 0);
        iStruct = 0;
        for i = 1:size(inCell, 1)
            if (inCell{i, 1}(1) == '#')
                iStruct = iStruct + 1;
                outStruct = [outStruct, struct()];
            else
                if (inCell{i, 3}(1) == '#')
                    outStruct{iStruct}.(inCell{i, 1}) = origCell{i, 4};
                else
                    outStruct{iStruct}.(inCell{i, 1}) = inCell{i, 2};
                end
            end
        end
    end

    function cellStr = num2cellStr(numVector)
        cellStr = cell(size(numVector));
        for iX = 1:numel(numVector)
            cellStr{iX} = num2str(numVector(iX));
        end
    end

    function SaveFigureToFile(hX, fileName)
        % Save figure from axes with handle hX to file
        fig2save = figure('Visible', 'off');
        newAxes = copyobj(hX, fig2save);
        set(newAxes, 'Units', 'normalized', 'Position', [0.13, 0.11, 0.775, 0.815]);
        saveas(fig2save, fileName, 'fig');
        close(fig2save);
        % Change visibility to 'on'
        map = load(fileName, '-mat');
        names = fieldnames(map);
        for j = 1:numel(names)
            map.(names{j}).properties.Visible = 'on';
        end
        save(fileName, '-struct', 'map');
    end
end