classdef TimeParamsCl
    properties (Access = public)
        t;              % Time vector
        dt;             % Time step
        nIntervals;     % Number of intervals
        nDays;          % Number of days
        intervalsPerDay;% Intervals per day
    end
%     properties (Dependent)
%         
%     end
    
    methods (Access = public)
        function self = TimeParamsCl(varargin)
            switch numel(varargin)
                case 1
                    if isstruct(varargin{1})
                        % Received structure with properties
                        self.nIntervals = varargin{1}.nIntervals;
                        self.nDays = varargin{1}.nDays;
                        self.t = varargin{1}.t;
                        self.dt = varargin{1}.dt;
                        self.intervalsPerDay = varargin{1}.intervalsPerDay;
                    else
                        % Received time vector
                        self.t = varargin{1};
                        self.dt = self.t(2) - self.t(1);
                        self = self.updNumIntervals();
                        self = self.updIntervalsPerDay();
                        self = self.updNumDays();
                    end
                case 2
                    % Received num days and time step
                    self.nDays = varargin{1};
                    self.dt = varargin{2};
                    self.t = 0:self.dt:self.nDays;
                    self = self.updIntervalsPerDay();
                    self = self.updNumIntervals();
                otherwise
                    % Nothing so far
                    error('Not implemented');
            end
        end
        
        function self = setT(self, newT)
            self.t = newT;
            self.dt = self.t(2) - self.t(1);
            self = self.updNumIntervals();
            self = self.updIntervalsPerDay();
            self = self.updNumDays();
        end
        
        function self = setNDays(self, newDays)
            self.nDays = newDays;
            self.t = self.t(self.t <= newDays);
            self = self.updNumIntervals();
        end
        
        function self = setNIntervals(self, newIntervals)
            self.nIntervals = newIntervals;
            self.t = self.t(1:self.nIntervals);
            self.nDays = updNumDays();
        end
    end
    
    methods (Access = private)
        function self = updIntervalsPerDay(self)
            self.intervalsPerDay = 1 / self.dt;
        end
            
        function self = updNumDays(self)
            self.nDays = self.nIntervals / self.intervalsPerDay;
        end
        
        function self = updNumIntervals(self)
            self.nIntervals = numel(self.t);
        end
    end
end