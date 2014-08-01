function result = RunOrchestraInterface(chemistryFilePath, varDefinitionTable, ioVariableList, inputVar, inputVarVal)
    % Run ORCHESTRA model using MATLAB-Java interface for ORCHESTRA.
    % Inputs:
    %   * chemistryFilePath - a path to file defining a chemistry for the
    %   model;
    %   * varDefinitionTable - a N-by-2 cell array defining names of N
    %   model variables (as char row-arrays) in the first column and their
    %   default values (numeric values);
    %   * ioVariableList - a M-element column-cell-vector defining a list
    %   of IO variables for the model;
    %   * inputVar - name or index of an IO variable in ioVariableList
    %   vector to loop over (if index, it must be in a range from 1 to M);
    %   * inputVarVal - 1-D vector of K values for the IO variable with
    %   index inputVarIdx.
    % Output:
    %   result - K-by-M array of values for M IO variables corresponding to
    %   K different values of an IO variable with index inputVarIdx.

    if ischar(inputVar)
        inputVarIdx = -1;
        for i = 1:numel(ioVariableList)
            if (strcmp(ioVariableList{i}, inputVar))
                inputVarIdx = i;
            end
        end
        if (inputVarIdx < 0)
            error('%s - no such IO variable!', inputVar);
        end
    elseif isnumeric(inputVar)
        inputVarIdx = inputVar;
    else
        error('Wrong format for setting IO variable!');
    end
    
    % First we have to change current directory to the directory containing
    % OrchestraInterface.jar file and other *.txt files with configuration
    % (we change back at the end of this cript)
    [stack, ~] = dbstack('-completenames');
    [orchestraDir, ~, ~] = fileparts(stack(1).file);
    workingDir = cd(orchestraDir);
    addpath(workingDir);
    
    try     % If error occures we change back to the original directory
    
        % Specify the path to *.inp file which defined our model's chemistry
        % (relative to this Matlab script, or absolute)
        chemistryFileName = [workingDir chemistryFilePath];
        
        % Next we specify the path to a compiled JAR file with the ORCHESTRA
        % interface. And import corresponding libraries
        javaclasspath([orchestraDir '/OrchestraInterface.jar']);
        import OrchestraInterface.*;
        
        % We initialize a list of variables
        numVars = size(varDefinitionTable, 1);
        ioIndices = ~ismember(ioVariableList, varDefinitionTable(:, 1));
        ioVariableListToNode = ioVariableList(ioIndices);
        numIoVars = numel(ioVariableList);
        numIoVarsToNode = numel(ioVariableListToNode);
        
        variableList = NodeVariableInfoList(numVars + numIoVarsToNode);
%         variableList = NodeVariableInfoList(numVars);
        
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

        % Initialize results
        result = zeros(numIoVars, numel(inputVarVal));        
        
        for idx = 1:numel(inputVarVal)
            % Run calculation method. it's inputs are:
            % - vector of indices of IO variables to be set (as in variable list);
            % - vector of values to be set for such variables;
            result(:, idx) = orchestraInstance.Calculate(inputVarIdx, inputVarVal(idx));
        end
        
    catch err
        
        cd(workingDir);
        rethrow(err);
        
    end
        
    % Go back to original working directory
    cd(workingDir);
    
    result = result';
end