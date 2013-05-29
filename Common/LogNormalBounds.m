function confBounds = LogNormalBounds(mu, sigma, nSigma)
%
%   This function returns an array of two elements which are the lower and upper bounds of
%   confidence interval for a log-normally distributed variable with given parameters mu and
%   sigma.
%       nSigma - defines the width of confidence interval like for normal distribution:
%           1   - 68.3%
%           2   - 95.4%
%           3   - 99.7%
%           etc.
%

    confBoundsNorm = [mu - nSigma * sigma, mu + nSigma * sigma];
    confBounds = exp(confBoundsNorm);
end