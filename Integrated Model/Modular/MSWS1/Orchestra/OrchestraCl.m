classdef OrchestraCl < ComponentCl
    properties (Access = private)
        isInstantiated = false;
        Ori;            % Orchestra instance(s)
        variableList;   % Orchestra variables
        Comp;
    end
    
    methods (Access = public)
        % Constructor is private to prevent multiple instantiations
        function self = OrchestraCl(Comp)
            % Specify the path to a compiled JAR file with the ORCHESTRA interface.
            % And import corresponding libraries
            ROOT_FOLDER = [cd '/../MSWS1/Orchestra'];
            ORCHESTRA_JAR = [ROOT_FOLDER '/Interface/OrchestraInterface.jar'];
            javaclasspath('-v1');
            jp = javaclasspath;
            jpFull = cell(size(jp));
            for i = 1:numel(jp)
                jpFull{i} = which(jp{i});
            end
            if ~ismember(which(ORCHESTRA_JAR), jpFull)
                warning(['Adding Orchestra to Matlab path. This may cause your program to misbehave. ' ...
                    'In such case stry restarting your program without rebooting Matlab.']);
                javaaddpath(ORCHESTRA_JAR);
            end
            import OrchestraInterface.*;
            
            % Initialize a list of Master species and I/O variables
            % Give each variable a name, default value, indicate if it is a static variable and indicate where this variable was defined.
            self.variableList = NodeVariableInfoList(length(Comp.all));
            ioVariableList = cell(1, length(Comp.all));
            for i = 1:length(Comp.all)
                self.variableList.SetVariable(i, Comp.all(i), Comp.alli(i), false, 'defined by example');
                ioVariableList(i) = Comp.all(i);
            end

            self.Ori = OrchestraModule([ROOT_FOLDER '/Bioreactor/chemistry.inp'], ...
                self.variableList, ioVariableList);
            self.Comp = Comp;
            self.isInstantiated = true;
        end
        
        function concentration = Calculate(self, varargin)
            % Initializes ORCHESTRA with total/derived concentrations in Bioreactor/chemistry.inp, returns JAVA module
            if nargin == 1
                concentration = self.Ori.Calculate(1:length(self.Comp.masteri), [self.Comp.masteri]);
            else
                concentration = self.Ori.Calculate(varargin{:});
            end
        end
    end
end