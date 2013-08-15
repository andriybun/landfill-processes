function ORI = initialize_ORI(Comp,x)
% Specify the path to a compiled JAR file with the ORCHESTRA interface.
% And import corresponding libraries
javaclasspath([cd '/../Orchestra/Interface/OrchestraInterface.jar']);
import OrchestraInterface.*;

% Initialize a list of Master species and I/O variables
% Give each variable a name, default value, indicate if it is a static variable
% and indicate where this variable was defined.
    variableList = NodeVariableInfoList(length(Comp.all));
    for i = 1:length(Comp.all)
        variableList.SetVariable(i, Comp.all(i), Comp.alli(i), false, 'defined by example');
        ioVariableList(i) = Comp.all(i);
    end

% Initialize ORCHESTRA object
% Initialize H2CO3.tot and H+.tot with initial pH2CO3 and pH
    if x ~= 0
        k1 = find(strcmp('H2CO3.logact',Comp.all));
        k2 = find(strcmp('pH',Comp.all));
        ORI = OrchestraModule([cd '/../Orchestra/Initialize/chemistry.inp'], variableList, ioVariableList);
        ORI = ORI.Calculate([k1 k2], [Comp.alli([k1 k2])]);
        clc
% Initialize ORCHESTRA with total amounts of all master species
    else
        ORI = OrchestraModule([cd '/../Orchestra/Bioreactor/chemistry.inp'], variableList, ioVariableList);
        
        activate_activity = ORI.Calculate(1:12, [Comp.masteri]);
        clc
    end
end