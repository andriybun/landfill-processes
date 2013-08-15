classdef progressbar
    properties
        maxVal;
        pbHandle;
    end
    
    methods
        function self = progressbar(maxVal, pbName)
            if nargin == 0
                maxVal = 1;
            end
            self.maxVal = maxVal;
            if nargin < 2
                pbName = '';
            end
            self.pbHandle = waitbar(0, sprintf('%3.2f%% done', 0), 'Name', pbName);
        end
        
        function self = update(self, newProg)
            newProg = newProg / self.maxVal;
            waitbar(newProg, self.pbHandle, sprintf('%3.2f%% done', newProg * 100));
        end
        
        function delete(self)
            close(self.pbHandle);
        end
    end
end