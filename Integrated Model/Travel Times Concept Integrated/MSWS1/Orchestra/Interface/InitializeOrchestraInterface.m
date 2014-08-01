function [orchestraInstance, loopVarsOut] = InitializeOrchestraInterface(chemistryFilePath, varDefinitionTable, ioVariableList, loopVars)
    % Run ORCHESTRA model using MATLAB-Java interface for ORCHESTRA.
    % Inputs:
    %   * chemistryFilePath - a path to file defining a chemistry for the
    %   model;
    %   * varDefinitionTable - a N-by-2 cell array defining names of N
    %   model variables (as char row-arrays) in the first column and their
    %   default values (numeric values);
    %   * ioVariableList - a M-element column-cell-vector defining a list
    %   of IO variables for the model;
    %   * loopVars - cell array of names of IO variables loop over.
    % Output:
    %   * orchestraInstance - an instance of ORCHESTRA interface object
    %   which is initialized and will be used later for computations;
    %   * loopVarsOut - array of structures with fields 'name', 'index'
    %   storing names of loop variables and their indices to be used for
    %   setting values for calculations.

    loopVarsOut = struct('name', {}, 'index', {});
    loopVarsOut = repmat(loopVarsOut, [1, numel(loopVars)]);
    
    for varIdx = 1:numel(loopVars)
        loopVarsOut(varIdx).name = loopVars{varIdx};
        loopVarsOut(varIdx).index = -1;
        for i = 1:numel(ioVariableList)
            if (strcmp(ioVariableList{i}, loopVarsOut(varIdx).name))
                loopVarsOut(varIdx).index = i;
            end
        end
        if (loopVarsOut(varIdx).index < 0)
            error('%s - no such IO variable!', loopVarsOut(varIdx).name);
        end
    end
    
    % First we have to change current directory to the directory containing
    % OrchestraInterface.jar file and other *.txt files with configuration
    % (we change back at the end of this cript)
    [stack, ~] = dbstack('-completenames');
    [orchestraDir, ~, ~] = fileparts(stack(1).file);
%     addpath(orchestraDir);
    
    try     % If error occures we change back to the original directory
    
        % Specify the path to *.inp file which defined our model's chemistry
        % (relative to this Matlab script, or absolute)
        chemistryFileName = [cd chemistryFilePath];
        
        % Next we specify the path to a compiled JAR file with the ORCHESTRA
        % interface. And import corresponding libraries
        javaclasspath([orchestraDir '/OrchestraInterface.jar']);
        import OrchestraInterface.*;
        
        % We initialize a list of variables
        numVars = size(varDefinitionTable, 1);
        ioIndices = ~ismember(ioVariableList, varDefinitionTable(:, 1));
        ioVariableListToNode = ioVariableList(ioIndices);
        numIoVarsToNode = numel(ioVariableListToNode);
        
        variableList = NodeVariableInfoList(numVars + numIoVarsToNode);
        
        % Then we define the set of variables that is stored in each Node of
        % this type. We give each variable a name, default value, indicate if
        % it is a static variable (all nodes share same variable) and indicate
        % where this variable was defined.
        for idx = 1:numVars
            variableList.SetVariable(idx, varDefinitionTable{idx, 1}, varDefinitionTable{idx, 2}, false, 'defined by example');
        end
        % We also add to the node IO variables
        for idx = 1:numIoVarsToNode
            variableList.SetVariable(numVars + idx, ioVariableListToNode{idx}, 0, false, 'defined by example');
        end
        
        % Initialize ORCHESTRA object
        orchestraInstance = OrchestraModule(chemistryFileName, variableList, ioVariableList);

    catch err
        
        rethrow(err);
        
    end
        
end