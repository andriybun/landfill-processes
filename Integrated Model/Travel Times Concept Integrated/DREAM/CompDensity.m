function [p,log_p] = CompDensity(x,MCMCPar,Measurement,ModelName,Extra,option)
% This function computes the density of each x value

p = []; log_p = [];

% Loop over the individual parameter combinations of x
for ii = 1:size(x,1),
    % Call model to generate simulated data
    evalstr = ['ModPred = ',ModelName,'(x(ii,:),Extra);']; eval(evalstr);
    
    if option == 1, % Model directly computes posterior density
        p(ii,1:2) = [ModPred ii]; log_p(ii,1) = log(p(ii,1));
    end;

    if option == 2, % Model computes output simulation
        Err = (Measurement.MeasData(:)-ModPred(:));
        % Compute the number of measurement data
        N = size(Measurement.MeasData,1);
        % Derive the log likelihood
        log_p(ii,1) = N.*log(MCMCPar.Wb./Measurement.Sigma) - MCMCPar.Cb.*(sum((abs(Err./Measurement.Sigma)).^(2/(1+MCMCPar.Gamma))));        
        % And retain in memory
        p(ii,1:2) = [log_p(ii,1) ii]; 
    end;

    if option == 3, % Model computes output simulation
        Err = (Measurement.MeasData(:)-ModPred(:));
        % Derive the sum of squared error
        SSR = sum(abs(Err).^(2/(1+MCMCPar.Gamma)));
        % And retain in memory
        p(ii,1:2) = [-SSR ii]; log_p(ii,1) = -0.5 * SSR;
    end;

    if option == 4, % Model directly computes log posterior density
        p(ii,1:2) = [ModPred ii]; log_p(ii,1) = p(ii,1);
    end;

    if option == 5, % Similar as 3, but now weights with the Measurement Sigma
        % Defime the error
        Err = (Measurement.MeasData(:)-ModPred(:));
        % Derive the sum of squared error
        SSR = sum(abs(Err).^(2/(1+MCMCPar.Gamma)));
        % And retain in memory
        p(ii,1:2) = [-SSR ii]; log_p(ii,1) = -0.5 * SSR;
    end;

    if option == 8, % Generalized log likelihood (GL)
        % Extract statistical model parameters
        par = Extra.fpar;               % fixed parameters
        par(Extra.idx_vpar) = x(ii,:);  % variable parameters
        par = par';                     % make it a column vector
        statpar = par(end-10:end);
        % Compute the log-likelihood
        log_p(ii,1) = GL('est',statpar,ModPred,Measurement.MeasData);
        % And retain in memory
        p(ii,1:2) = [log_p(ii,1) ii];
    end;

end;