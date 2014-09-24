function [p,log_p] = CompDensity_par(x,MCMCPar,Measurement,ModelName,Extra,option)
% This function computes the density of each x value

% To paralellize function
numWorkers = Extra.numWorkers;
[ptmp, log_p] = paralellize(@f, x, numWorkers, MCMCPar, Measurement, ModelName, Extra, option);
p = [ptmp' (1:length(ptmp))']; log_p = log_p';
clc;
return
end

function [ptmp, log_p] = f(x_in, MCMCPar, Measurement, ModelName, Extra, option)

% Re-initialize ORCHESTRA
Extra.ORI = initialize_ORI(Extra.Comp,0);

% Call model to generate simulated data
evalstr = ['ModPred = ',ModelName,'(x_in,Extra);']; eval(evalstr);

if option == 1, % Model directly computes posterior density
    ptmp = ModPred; log_p = log(ptmp); 
end

if option == 2, % Model computes output simulation
    Err = (Measurement.MeasData(:)-ModPred(:));
    % Compute the number of measurement data
    N = size(Measurement.MeasData,1);
    % Derive the log likelihood
    log_p(ii,1) = N.*log(MCMCPar.Wb./Measurement.Sigma) - MCMCPar.Cb.*(sum((abs(Err./Measurement.Sigma)).^(2/(1+MCMCPar.Gamma))));
    % And retain in memory
    ptmp = log_p;
end;

if option == 3, % Model computes output simulation
    Err = (Measurement.MeasData(:)-ModPred(:));
    % Derive the sum of squared error
    SSR = sum(abs(Err).^(2/(1+MCMCPar.Gamma)));
    % And retain in memory
    ptmp = -SSR; log_p = -0.5 * SSR;
end;

if option == 4, % Model directly computes log posterior density
    ptmp = ModPred; log_p = ptmp; 
end;

if option == 5, % Similar as 3, but now weights with the Measurement Sigma
    % Defime the error
    Err = (Measurement.MeasData(:)-ModPred(:));
    % Derive the sum of squared error
    SSR = sum(abs(Err).^(2/(1+MCMCPar.Gamma)));
    % And retain in memory
    ptmp = -SSR; log_p = -0.5 * SSR;  
end;

if option == 8, % Generalized log likelihood (GL)
    % Extract statistical model parameters
    par = Extra.fpar;               % fixed parameters
    par(Extra.idx_vpar) = x(ii,:);  % variable parameters
    par = par';                     % make it a column vector
    statpar = par(end-10:end);
    % Compute the log-likelihood
    log_p = GL('est',statpar,ModPred,Measurement.MeasData);
    % And retain in memory
    ptmp = log_p;
end;
end