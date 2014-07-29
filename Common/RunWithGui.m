function output = GenerateGui(input, func)

    inputCell = struct2cell(input);
    CreateGui(inputCell, func);
    
    output = 3;
    
    return
    
    function CreateGui(inputCell, func)
        % Window properties
        WINDOW_WIDTH = 400;
        WINDOW_HEIGHT = 500;
        % Progress bar properties
        PRBAR_HEIGHT = 20;
        PRBAR_MARGIN = 10;
        % Table properties
        MAX_TABLE_HEIGHT = 200;
        COLUMN_WIDTHS_REL = [70, 180, 50];
        ROW_HEIGHT = 19;
        BOUND_PIX = 2;
        % Button propeties
        BUTTON_WIDTH = 70;
        BUTTON_HEIGHT = 25;
        % Pop-up properties
        LIST_HEIGHT = 20;
        
        % Create and then hide the GUI as it is being constructed.
        guiWindow = figure('Visible', 'off', 'Position', [200, 200, WINDOW_WIDTH, WINDOW_HEIGHT]);
        set(guiWindow, 'Resize', 'off');
        
        % Progress bar
        prBar = progressbar(guiWindow, ...
            [PRBAR_MARGIN, PRBAR_MARGIN, WINDOW_WIDTH - 2 * PRBAR_MARGIN, PRBAR_HEIGHT], [0 100]);
        
        % Construct components
        % Button to run model
        btnVertPos = WINDOW_HEIGHT - 10 - BUTTON_HEIGHT;
        
        % Table
        tableWidth = WINDOW_WIDTH - 2 * PRBAR_MARGIN;
        columnWidths = round(tableWidth / sum(COLUMN_WIDTHS_REL) * COLUMN_WIDTHS_REL);
        columnWidths = columnWidths - [BOUND_PIX, sum(columnWidths) - tableWidth - 1, BOUND_PIX];
        tableHeight = min(MAX_TABLE_HEIGHT, ...
            (size(inputCell, 1) + 1) * ROW_HEIGHT + BOUND_PIX * 2 - 1);
        tableVertPos = btnVertPos - 10 - tableHeight;
        hTable = uitable('Position', ...
            [10, tableVertPos, tableWidth, tableHeight]);
        set(hTable, 'ColumnName', {'Parameter', 'Value', 'Type'});
        set(hTable, 'ColumnEditable', [false, true, false]);
        set(hTable, 'ColumnWidth', num2cell(columnWidths));
        set(hTable, 'RowName', []);
        set(hTable, 'Data', inputCell);
        
        % Pop-up list
        popListVertPos = tableVertPos - 10 - LIST_HEIGHT;
        popListWidth = (WINDOW_WIDTH - 3 * PRBAR_MARGIN) / 2;
        hPopListX = uicontrol('Style', 'popupmenu', 'String', 'x-axis data', 'Enable', 'off', ...
            'Position', [PRBAR_MARGIN, popListVertPos, popListWidth, LIST_HEIGHT]);
        popListYhorOffset = 2 * PRBAR_MARGIN + popListWidth;
        hPopListY = uicontrol('Style', 'popupmenu', 'String', 'y-axis data', 'Enable', 'off', ...
            'Position', [popListYhorOffset, popListVertPos, popListWidth, LIST_HEIGHT]);
        
        % Plot
        plotHeight = WINDOW_HEIGHT - popListVertPos - PRBAR_HEIGHT - 10;
        plotWidth = WINDOW_WIDTH - 6 * PRBAR_MARGIN;
        plotHorPos = WINDOW_WIDTH - PRBAR_MARGIN - plotWidth;
        plotVertPos = popListVertPos - 2 * 10 - plotHeight;
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
        set(hPopListX, 'Callback', {@PlotGraph, uiElements});
        set(hPopListY, 'Callback', {@PlotGraph, uiElements});
        
        hBtnStart = uicontrol('Style', 'pushbutton', ...
            'String', 'Start', ...
            'Tag', 'hBtnStart', ...
            'Position', [10, btnVertPos, BUTTON_WIDTH, BUTTON_HEIGHT], ...
            'Callback', {@FuncWrap, {func, inputCell, uiElements}});
        
%         ha = axes('Units','Pixels','Position',[50,60,200,185]);
%         align([hsurf,hmesh,hcontour,htext,hpopup],'Center','None');
        
        %Make the GUI visible.
        set(guiWindow, 'Visible', 'on');
    end
    
    function output = FuncWrap(hObject, eventdata, args)
        output = args{1}(args{2}, args{3}.prBar);
        set(args{3}.hPopListX, 'String', fieldnames(output));
        set(args{3}.hPopListX, 'Enable', 'on');
        set(args{3}.hPopListY, 'String', fieldnames(output));
        set(args{3}.hPopListY, 'Enable', 'on');
        set(hObject, 'UserData', output);
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

    function outCell = struct2cell(inStruct)
        %% TODO: do something with non-single value parameters
        fieldNameVec = fieldnames(inStruct);
        nFields = numel(fieldNameVec);
        outCell = cell(nFields, 3);
        for i = 1:nFields
            fieldName = fieldNameVec{i};
            fieldVal = inStruct.(fieldName);
            outCell(i, :) = {fieldName, fieldVal, class(fieldVal)};
        end
    end
end