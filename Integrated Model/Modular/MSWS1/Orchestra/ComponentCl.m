classdef ComponentCl
    methods (Abstract)
        concentration = Calculate(self);
    end
end