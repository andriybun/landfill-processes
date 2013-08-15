function [chemistryFilePath, varDefinitionTable, ioVariableList] = ParseOrchestraModel(modelDir)
    % Parse ORCHESTRA model files and prepare inputs for running model from
    % MATLAB

    %% Parse composer.inp with main file names
    fileNames = ParseComposer([modelDir 'composer.inp']);
    
    %% Parse colum file with basic varialbes
    iniVarDefinitionTable = ParseColumn([modelDir fileNames.columndat]);

    chemistryFilePath = ['/' modelDir fileNames.chemistry];

    %% Parse file with variables and their default values
    [varList, varVals] = ReadVariables([modelDir fileNames.input]);
    varDefinitionTable = [varList varVals];
    varDefinitionTable = vertcat(iniVarDefinitionTable, varDefinitionTable);
    
    %% Parse file with the list of IO variables
    [ioVariableList, ~] = ReadVariables([modelDir fileNames.outputformat]);
    
    return

    %% Functions for parsing model parameters
    function files = ParseComposer(fileName)
        composerFile = textread(fileName, '%s', 'delimiter', '\n');    
        
        varsList = { };
        
        files = struct();
        
        for i = 1:numel(composerFile)
            line = composerFile{i};
            
            % Skip lines that are: empty, commments
            if (isempty(line) || (strcmp(line(1:2), '//')))
                continue
            end
            
            % Split strings to a separate words
            tk = TokenizeString(line);
            
            if (strcmp(tk{1}, 'HFile:') && strcmp(tk{2}, 'outputformat'))
                files.outputformat = tk{end};
            end
            
            if (strcmp(tk{1}, 'CHEMFile:') && strcmp(tk{2}, 'Chemistry'))
                files.chemistry = tk{end};
            end
            
            if (strcmp(tk{1}, 'HFile:') && strcmp(tk{2}, 'Column'))
                files.columndat = tk{end};
            end
            
            if (strcmp(tk{1}, 'File:') && strcmp(tk{2}, 'Input'))
                files.input = tk{end};
            end
        end
    end
    
    function varsList = ParseColumn(fileName)
        columnFile = textread(fileName, '%s', 'delimiter', '\n');    
        
        varsList = { };
        
        for i = 1:numel(columnFile)
            line = columnFile{i};
            
            % Skip lines that are: empty, commments
            if (isempty(line) || (strcmp(line(1:2), '//')))
                continue
            end
            
            % Split strings to a separate words
            tk = TokenizeString(line);
            
            if (strcmp(tk{1}, '@Var:'))
                varsList = vertcat(varsList, [tk(2), str2double(tk(3))]);
            end
        end
    end

    function [varList, varVals] = ReadVariables(fileName)
        varList = { };
        varVals = { };
        
        inpFile = textread(fileName, '%s', 'delimiter', '\n');
        
        for i = 1:numel(inpFile)
            line = inpFile{i};
            
            % Skip lines that are: empty, commments
            if (isempty(line) || (strcmp(line(1:2), '//')))
                continue
            end
            
            % Split strings to a separate words
            tk = TokenizeString(line);
            
            if (strcmp(tk{1}, 'Var:'))
                varList = tk(2:end)';
            end
            
            if (strcmp(tk{1}, 'Data:'))
                varVals = num2cell(str2double(tk(2:end)'));
            end
        end
    end

    function tokens = TokenizeString(inStr)
        [tokens, rem] = strtok(inStr);
        while ~isempty(rem)
            [tok, rem] = strtok(rem);
            if ~isempty(tok)
                tokens = [tokens, {tok}];
            end
        end
    end

end