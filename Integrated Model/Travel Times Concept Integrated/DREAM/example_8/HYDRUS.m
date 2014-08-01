function [log_likelihood] = HYDRUS(x,Extra)
% Runs the HYDRUS model and computes the log-likelihood

% Back-transform parameters
par = x; par(3:5) = 10.^par(3:5);

try

    % Run HYDRUS-1D
	sims = runH1D(par,Extra);
	
	% Filter sims for measurement dates
	ind = zeros(size(Extra.hoy));
	for i = 1:size(Extra.hoy,1)
		ind(i) = find(sims.hoy == Extra.hoy(i));
	end

	% Compute number of measurements
	N = length(sims.water(ind));
	
	% Compute residuals
	epsilon = Extra.water - sims.water(ind);
	
	% Compute log-likelihood
	log_likelihood = - (N/2)*log(sum(epsilon.^2));
	
catch	
	% If HYDRUS-1D did not converge, set the log-likelihood to some arbitrary low value
	log_likelihood = -9999;
end

log_likelihood