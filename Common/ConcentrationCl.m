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
        
        function self = SetVolume(self, v, varargin)
            if (nargin == 2)
                self.v = v;
            elseif (nargin > 2)
                self.v(varargin{:}) = v(varargin{:});
            end
        end
        
        function self = SetConcentration(self, c, varargin)
            if (nargin == 2)
                self.m = self.v .* c;
            elseif (nargin > 2)
                self.m(varargin{:}) = self.v(varargin{:}) .* c;
            end
        end
        
        function self = AddVolume(self, dv, varargin)
            if (nargin == 2)
                self.v = self.v + dv;
            elseif (nargin > 2)
                self.v(varargin{:}) = self.v(varargin{:}) + dv;
            end
        end
        
        function self = AddSolute(self, dv, dc)
            self = self.AddVolume(dv);
            dm = dv .* dc;
            self.m = self.m + dm;
        end
        
        function v = GetVolume(self, varargin)
            if (nargin == 1)
                v = self.v;
            elseif (nargin > 1)
                v = self.v(varargin{:});
            end
        end
        
        function c = GetConcentration(self, varargin)
            if (nargin == 1)
                c = self.m ./ self.v;
                isZero = (self.v == 0);
            elseif (nargin > 1)
                c = self.m(varargin{:}) ./ self.v(varargin{:});
                isZero = (self.v(varargin{:}) == 0);
            end
            c(isZero) = 0;
        end
        
        function m = GetMass(self, varargin)
            if (nargin == 1)
                m = self.m;
            elseif (nargin > 1)
                m = self.m(varargin{:});
            end
        end
    end
    
    methods (Access = private)
        
    end
    
end