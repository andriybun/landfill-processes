function output = RunWithGui(appName, input, func)
% This function is used for running generic function using input parameters defined. It is also able
% to plot the result.
% Inputs:
%   appName - string. The GUI window will have the name as specified by this parameter.
%   input - a structure with parameters, that will be passed to the function. Note, composite types
%           of fields (vectors, matrices, cells, sub-structures) may not work properly. In any case,
%           it will be not possible to edit values of them.
%   func - handle of a function to run. The function may accept the inputs passed as 'input' 
%           structure and a progress bar [0, 100]. Its values can be set using 'pvalue' property.
%           The function must return a structure with results. Its fields then can be plotted using
%           GUI.
% This function generates GUI window that allows change values for properties contained in input
% structure.
% After the function defined by handle finishes, a user can select its output data to be plot. For
% this the function should return a structure with vectors for each data entity.

    inputCell = struct2cellX(input);
    CreateGui(inputCell, func);
    
    output = 3;
    
    return
    
    function CreateGui(inputCell, func)
        % Window properties
        WINDOW_WIDTH = 400;
        WINDOW_HEIGHT = 600;
        % Progress bar properties
        PRBAR_HEIGHT = 20;
        PRBAR_MARGIN = 10;
        % Table properties
        MAX_TABLE_HEIGHT = 200;
        COLUMN_WIDTHS_REL = [70, 160, 70];
        ROW_HEIGHT = 19;
        BOUND_PIX = 2;
        % Button propeties
        BUTTON_WIDTH = 70;
        BUTTON_HEIGHT = 25;
        % Pop-up properties
        LIST_HEIGHT = 20;
        
        % Create and then hide the GUI as it is being constructed.
        guiWindow = figure('Visible', 'off', 'Position', [200, 200, WINDOW_WIDTH, WINDOW_HEIGHT]);
        set(guiWindow, 'Name', appName, 'NumberTitle', 'off', 'Tag', 'guiWindow');
        set(guiWindow, 'Resize', 'off');
        
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
        
        % Pop-up list
        popListVertPos = tableVertPos - PRBAR_MARGIN - LIST_HEIGHT;
        popListWidth = (WINDOW_WIDTH - 3 * PRBAR_MARGIN) / 2;
        hPopListX = uicontrol('Style', 'popupmenu', 'String', 'x-axis data', 'Enable', 'off', ...
            'Position', [PRBAR_MARGIN, popListVertPos, popListWidth, LIST_HEIGHT]);
        popListYhorOffset = 2 * PRBAR_MARGIN + popListWidth;
        hPopListY = uicontrol('Style', 'popupmenu', 'String', 'y-axis data', 'Enable', 'off', ...
            'Position', [popListYhorOffset, popListVertPos, popListWidth, LIST_HEIGHT]);
        
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
        
        uiElements = struct();
        uiElements.prBar = prBar;
        uiElements.hPopListX = hPopListX;
        uiElements.hPopListY = hPopListY;
        uiElements.hPlotAxes = hPlotAxes;
        uiElements.hTable = hTable;
        set(hPopListX, 'Callback', {@PlotGraph, uiElements});
        set(hPopListY, 'Callback', {@PlotGraph, uiElements});
        
        hBtnStart = uicontrol('Style', 'pushbutton', ...
            'String', 'Start', ...
            'Tag', 'hBtnStart', ...
            'Position', [10, btnVertPos, BUTTON_WIDTH, BUTTON_HEIGHT], ...
            'Callback', {@FuncWrap, {func, inputCell, uiElements}});
        
        % Make the GUI visible.
        set(guiWindow, 'Visible', 'on');
    end
    
    function FuncWrap(hObject, eventdata, args)
        func = args{1};
        origInputs = args{2};
        uiElements = args{3};
        output = func(cell2structX(get(uiElements.hTable, 'data'), origInputs), uiElements.prBar);
        set(uiElements.hPopListX, 'String', fieldnames(output));
        set(uiElements.hPopListX, 'Enable', 'on');
        set(uiElements.hPopListY, 'String', fieldnames(output));
        set(uiElements.hPopListY, 'Enable', 'on');
        set(hObject, 'UserData', output);
        cla(uiElements.hPlotAxes);
        guidata(hObject);
    end

    function PlotGraph(hObject, eventdata, uiElements)
        output = get(findobj('Tag', 'hBtnStart'), 'UserData');
        xStrings = get(uiElements.hPopListX, 'String');
        xAxisData = xStrings{get(uiElements.hPopListX, 'Value')};
        yStrings = get(uiElements.hPopListY, 'String');
        yAxisData = yStrings{get(uiElements.hPopListY, 'Value')};
        plot(uiElements.hPlotAxes, output.(xAxisData), output.(yAxisData));
        xlabel(uiElements.hPlotAxes, xAxisData);
        ylabel(uiElements.hPlotAxes, yAxisData);
    end

    function outCell = struct2cellX(inStruct)
        fieldNameVec = fieldnames(inStruct);
        nFields = numel(fieldNameVec);
        outCell = cell(nFields, 4);
        for i = 1:nFields
            fieldName = fieldNameVec{i};
            fieldVal = inStruct.(fieldName);
            sz = size(fieldVal);
            if isstruct(fieldVal)
                valStr = sprintf('struct %d x %d', sz(1), sz(2));
                outCell(i, :) = {fieldName, valStr, '#Struct', fieldVal};
            else
                if isequal(sz, [1, 1]) || ischar(fieldVal)
                    outCell(i, 1:3) = {fieldName, fieldVal, class(fieldVal)};
                else
                    valStr = sprintf('vector %d x %d', sz(1), sz(2));
                    classStr = sprintf('#%s', class(fieldVal));
                    outCell(i, :) = {fieldName, valStr, classStr, fieldVal};
                end
            end
        end
    end

    function outStruct = cell2structX(inCell, origCell)
        outStruct = struct();
        for i = 1:size(inCell, 1)
            if (inCell{i, 3}(1) == '#')
                outStruct.(inCell{i, 1}) = origCell{i, 4};
            else
                outStruct.(inCell{i, 1}) = inCell{i, 2};
            end
        end
    end
end