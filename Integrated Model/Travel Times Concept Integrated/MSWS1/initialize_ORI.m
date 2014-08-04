function ORI = initialize_ORI(Comp, sw)
    % Specify the path to a compiled JAR file with the ORCHESTRA interface.
    % And import corresponding libraries
    ORCHESTRA_FOLDER = [cd '/../Orchestra/Interface/OrchestraInterface.jar'];
    javaclasspath('-v1');
    jp = javaclasspath;
    jpFull = cell(size(jp));
    for i = 1:numel(jp)
        jpFull{i} = which(jp{i});
    end
    if ~ismember(which(ORCHESTRA_FOLDER), jpFull)
        warning(['Adding Orchestra to Matlab path. This may cause your program to misbehave. ' ...
            'In such case stry restarting your program without rebooting Matlab.']);
        javaaddpath(ORCHESTRA_FOLDER);
    end
    
    % javaaddpath(ORCHESTRA_FOLDER);
    import OrchestraInterface.*;

    % Initialize a list of Master species and I/O variables
    % Give each variable a name, default value, indicate if it is a static variable and indicate where this variable was defined.
    variableList = NodeVariableInfoList(length(Comp.all));
    for i = 1:length(Comp.all)
        variableList.SetVariable(i, Comp.all(i), Comp.alli(i), false, 'defined by example');
        ioVariableList(i) = Comp.all(i);
    end
    
    switch sw
        % Initializes ORCHESTRA with total/derived concentrations in Bioreactor/chemistry.inp, returns JAVA module
        case 0
            ORI = OrchestraModule([cd '/./Orchestra/Bioreactor/chemistry.inp'], variableList, ioVariableList);
            activate_activity = ORI.Calculate(1:length(Comp.masteri), [Comp.masteri]); % activates activity calculation
            
            % Initializes ORCHESTRA with total/derived concentrations in Initialize/chemistry.inp, returns total & derived concentrations
        case 1
            ORI = OrchestraModule([cd '/./Orchestra/Initialize/chemistry.inp'], variableList, ioVariableList);
            ORI = ORI.Calculate(1:length(Comp.masteri), [Comp.masteri]);
            
            % Initializes ORCHESTRA with total/derived concentrations in Bioreactor/chemistry.inp, returns total & derived concentrations
        case 2
            ORI = OrchestraModule([cd '/./Orchestra/Bioreactor/chemistry.inp'], variableList, ioVariableList);
            ORI = ORI.Calculate(1:length(Comp.masteri), [Comp.masteri]);
    end
end