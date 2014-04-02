classdef LogNormalCl
    % Class for computing PDF, CDF, confidence intervals for Log-Normal distribution with delay
    properties (Access = public)
        mu;
        sigma;
        delay;
        bounds;
    end
    
    properties (Access = private)
        Const;
    end
    
    methods (Access = public)
        function self = LogNormalCl(mu, sigma, delay, Const)
            self.mu = mu;
            self.sigma = sigma;
            self.delay = delay;
            self.Const = Const;
            self.bounds = self.Bounds(self.mu, self.sigma, Const.NUM_SIGMAS) + self.delay;
        end
        
        function pdf = PdfDelayed(self, t)
            pdf = self.Pdf(t-self.delay, self.mu, self.sigma);
        end
    end
    
    methods (Static)
        function pdf = Pdf(t, mu, sigma)
            dt = t(end) - t(end-1);
            pdf = diff(logncdf([t, t(end) + dt], mu, sigma));
        end
        
        function confBounds = Bounds(mu, sigma, nSigma)
            %   This method returns an array of two elements which are the lower and upper bounds of
            %   confidence interval for a log-normally distributed variable with given parameters mu 
            %   and sigma.
            %       nSigma - defines the width of confidence interval like for normal distribution:
            %           1   - 68.3%
            %           2   - 95.4%
            %           3   - 99.7%
            %           etc.
            
            confBoundsNorm = [mu - nSigma * sigma, mu + nSigma * sigma];
            confBounds = exp(confBoundsNorm);
        end
    end
end
    