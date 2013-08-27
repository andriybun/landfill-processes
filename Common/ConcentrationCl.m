classdef ConcentrationCl
    
    properties (Access = public)
        % Number of elements
        nEl;
        % Number of solutes
        nSolutes;
    end
    
    properties (Access = private)
        % Number of dimensions
        nDims;
        % Volume of elements
        v;
        % Masses of solutes in each element (must have extra dimension for solutes, if nSolutes > 1)
        m;
    end
    
    methods (Access = public)
        function self = ConcentrationCl(vIni, cIni)
            cDim = size(cIni);
            self.nDims = numel(cDim);
            self.nSolutes = cDim(end);
            self.nEl = numel(vIni);
            self.v = vIni;
            self.m = self.op(vIni, cIni, @times);
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
                self.m(varargin{:}) = self.op(self.v(varargin{1:self.nDims-1}), c, @times);
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
            dm = self.op(dv, dc, @times);
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
                c = self.op(self.v(varargin{1:self.nDims-1}), self.m(varargin{:}), @self.divideInv);
                isZero = (self.v(varargin{1:self.nDims-1}) == 0);
            end
            nCopies = ones(1, self.nDims);
            nCopies(self.nDims) = numel(varargin{end});
            c(repmat(isZero, nCopies)) = 0;
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
        function res = op(self, ar1, ar2, fHandle)
            % Function is used to perform operation defined by function with handle fHandle on two
            % arrays, where second array has an extra dimension, while all other dimensions are the
            % same.
            arSize = size(ar2);
            nSolutesX = arSize(end);
            res = zeros(size(ar2));
            idx = cell(1, self.nDims);
            for iDim = 1:self.nDims-1
                idx{iDim} = ':';
            end
            for iSolute = 1:nSolutesX
                idx{self.nDims} = iSolute;
                res(idx{:}) = fHandle(ar1, ar2(idx{:}));
            end
        end
%     end
%     
%     methods (Static)
        function res = divideInv(self, v1, v2)
            res = v2 ./ v1;
        end
    end
    
end