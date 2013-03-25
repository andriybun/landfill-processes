classdef LogNormalParams
    properties
        OptParams;         % structure to which a file with best approximations of parameters will be loaded
        HydraulicParams;   % van Genuchten parameters used for simulations
        seLow;             % minimum value for precomputed effective saturation
        seHi;              % maximum value for precomputed effective saturation 
        logKsatRef;        % natural log of reference kSat
    end
    
    methods
        function self = LogNormalParams(file_name)
            if nargin == 0
                self.OptParams = load('OptParams_wt.mat');
            else
                self.OptParams = load(file_name);
            end
            if self.OptParams.saturation_effective_avg(1) > self.OptParams.saturation_effective_avg(end)
                self.seHi = self.OptParams.saturation_effective_avg(1);
                self.seLow = self.OptParams.saturation_effective_avg(end);
            else
                self.seLow = self.OptParams.saturation_effective_avg(1);
                self.seHi = self.OptParams.saturation_effective_avg(end);
            end
            self.HydraulicParams = self.OptParams.van_genuchten_params;
            self.OptParams = rmfield(self.OptParams, 'van_genuchten_params');
            self.logKsatRef = log(self.HydraulicParams.k_sat);
        end
        
        function [mu, sigma] = GetParams(self, kSat, se)
            mu = zeros(size(se));
            sigma = zeros(size(se));

            se = min(se, self.seHi);
            se = max(se, self.seLow);

            mu = self.logKsatRef - log(kSat) + interp1(self.OptParams.saturation_effective_avg, self.OptParams.mu, se, 'spline');
            sigma = interp1(self.OptParams.saturation_effective_avg, self.OptParams.sigma, se, 'spline');
        end
    end
end