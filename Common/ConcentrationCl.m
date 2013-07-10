classdef ConcentrationCl
    
    properties (Access = public)
        % Number of elements
        nEl;
        % Number of solutes
        nSolutes;
    end
    
    properties (Access = private)
        % Volume of elements
        v;
        % Masses of solutes in each element
        m;
    end
    
    methods (Access = public)
        function self = ConcentrationCl(vIni, cIni)
            self.nEl = numel(vIni);
            self.v = vIni;
            self.m = vIni .* cIni;
        end
        
        function self = SetVolume(self, v, i, j)
            if (nargin == 2)
                self.v = v;
            elseif (nargin == 3)
                self.v(i) = v(i);
            elseif (nargin == 4)
                self.v(i, j) = v(i, j);
            end
        end
        
        function self = SetConcentration(self, c, i, j)
            if (nargin == 2)
                self.m = self.v .* c;
            elseif (nargin == 3)
                self.m(i) = self.v(i) .* c;
            elseif (nargin == 4)
                self.m(i, j) = self.v(i, j) .* c;
            end
        end
        
        function self = AddVolume(self, dv, i, j)
            if (nargin == 2)
                self.v = self.v + dv;
            elseif (nargin == 3)
                self.v(i) = self.v(i) + dv;
            elseif (nargin == 4)
                self.v(i, j) = self.v(i, j) + dv;
            end
        end
        
        function self = AddSolute(self, dv, dc)
            self = self.AddVolume(dv);
            dm = dv .* dc;
            self.m = self.m + dm;
        end
        
        function v = GetVolume(self, i, j)
            if (nargin == 1)
                v = self.v;
            elseif (nargin == 2)
                v = self.v(i);
            elseif (nargin == 3)
                v = self.v(i, j);
            end
        end
        
        function c = GetConcentration(self, i, j)
            if (nargin == 1)
                c = self.m ./ self.v;
                isZero = (self.v == 0);
            elseif (nargin == 2)
                c = self.m(i) ./ self.v(i);
                isZero = (self.v(i) == 0);
            elseif (nargin == 3)
                c = self.m(i, j) ./ self.v(i, j);
                isZero = (self.v(i, j) == 0);
            end
            c(isZero) = 0;
        end
        
        function m = GetMass(self, i, j)
            if (nargin == 1)
                m = self.m;
            elseif (nargin == 2)
                m = self.m(i);
            elseif (nargin == 3)
                m = self.m(i, j);
            end
        end
    end
    
    methods (Access = private)
        
    end
    
end