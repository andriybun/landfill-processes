classdef OrchestraCompositeCl < ComponentCl
    properties (Access = private)
        children;
    end
    
    methods (Access = public)
        function self = OrchestraCompositeCl()
            self.children = {};
        end
        
        function self = add(self, newEl)
            if (~isa(newEl, 'ComponentCl'))
                error('Trying to add an object of wrong type');
            end
            self.children = cat(1, self.children, {newEl});
        end
        
        function result = Calculate(self)
            result = {};
            for i = 1:numel(self.children)
                result = cat(2, result, self.children{i}.Calculate());
            end
        end
    end
end