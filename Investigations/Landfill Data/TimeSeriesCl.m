classdef TimeSeriesCl
    
    properties
        headers;
        data;
        dateStart;
    end
    
    methods (Access = public)
        % Constructor
        function self = TimeSeriesCl(fileName)
            Raw = load(fileName);
            self.headers = Raw.headers;
            self.data = Raw.data;
            self.dateStart = self.GetTime(1);
        end
        
        %% Visualization
        function PrintHeaders(self)
            % Print headers and their indices
            [num2cell(1:numel(self.headers))', self.headers]
        end
        
        function Plot(self, iCols, tickMode)
            % Plot specified data series
            if nargin < 3
                tickMode = 'default';
            end
            switch tickMode
                case 'default'
                    plot(self.GetTime(), self.data(:, iCols));
                case 'date'
                    t = self.GetTime() - self.GetTime(1) + self.dateStart;
                    plot(t, self.data(:, iCols));
                    datetick('x', 'mmm yyyy', 'keepticks');
                otherwise
                    warning('Warning! Wrong tick mode. Plotting default ticks.');
                    plot(self.GetTime(), self.data(:, iCols));
            end
            legend(self.headers(iCols))
        end
        
        %% Data processing methods
        
        function self = ResetTime(self, tBase)
            % Reset time axis to start with tBase (if not provided, time starts with zero)
            if nargin < 2
                tBase = self.GetTime(1);
            end
            self.data(:, 1) = self.GetTime() - tBase;
        end
        
        function self = ChangeStartDate(self, newStDate)
            % Changes start date and if new date is after the original start date, the dataset is
            % cut off.
            self.data(:, 1) = self.GetTime() - self.GetTime(1) + self.dateStart - newStDate;
            self.data(self.GetTime() < 0, :) = [];
            self.dateStart = newStDate;
        end
        
        function t = GetTime(self, i)
            % Return time vector of dataset. If i is provided - returns only values indexed by it.
            if nargin < 2
                i = ':';
            end
            t = self.data(i, 1);
        end
        
        function self = InterpTime(self, argIn)
            % Interpolate values for a given time intervals/vector.
            % Note: this approach works only for cumulative datasets (levels also possible)!
            t = self.GetTime();
            if numel(argIn) == 1
                % If argIn is a single value - it's considered an interval
                tSamp = (t(1):argIn:t(end))';
            else
                % If argIn is a vector - it's considered to be a time sample for interpolation
                tSamp = argIn;
                self.dateStart = tSamp(1);
            end
            nCols = size(self.data, 2);
            dataSamp = nan(numel(tSamp), nCols-1);
            for iCol = 1:nCols-1
                iNotNan = ~isnan(self.data(:, iCol+1));
                dataSamp(:, iCol) = interp1(t(iNotNan), self.data(iNotNan, iCol+1), tSamp);
            end
            self.data = cat(2, tSamp, dataSamp);
        end
        
        function self = RemoveDrops(self, iCols)
            % Make cumulative data continuous, if for some reason values have been reset
            if nargin < 2
                iCols = 2:size(self.data, 2);
            end
            for iCol = iCols
                iNotNan = ~isnan(self.data(:, iCol));
                dataSample = self.data(iNotNan, iCol);
                diffSample = diff(dataSample);
                iDrops = find(diffSample < 0);
                for j = 1:numel(iDrops)
                    dataSample(iDrops(j)+1:end) = dataSample(iDrops(j)+1:end) - ...
                        diffSample(iDrops(j));
                end
                self.data(iNotNan, iCol) = dataSample;
            end
        end
        
        function self = DeleteCols(self, iCols)
            % Delete given columns from dataset
            self.headers(iCols) = [];
            self.data(:, iCols) = [];
        end
        
        function self = Join(self, other)
            % Join another dataset with the same time vector
            % First check if time vectors are equal
            if ~(isequal(self.GetTime(), other.GetTime()) && ...
                    (self.dateStart == other.dateStart))
                error('Sorry can''t join, time vectors and/or start dates are different');
            end
            self.headers = cat(1, self.headers, other.headers(2:end));
            self.data = cat(2, self.data, other.data(:, 2:end));
        end
    end
    
    
end
