classdef OrchestraCl < ComponentCl
    properties (Access = private)
        Ori;
    end
    
    methods (Access = public)
        % Constructor is private to prevent multiple instantiations
        function self = OrchestraCl(Ori)
            self.Ori = Ori;
            
        end
        
        function result = Calculate(self)
            result = self.Ori;
        end
    end
end