function [p,log_p] = CompDensity_par(x,MCMCPar,Measurement,ModelName,Extra,option)
% This function computes the density of each x value

[s,r] = size(x);
numWorkers = Extra.numWorkers;
[ptmp, log_p] = paralelize(@f, x, numWorkers, MCMCPar, Measurement, ModelName, Extra, option);
p = [ptmp' (1:length(ptmp))'];
log_p = log_p';
return

    function [ptmp, log_p] = f(x_in, MCMCPar, Measurement, ModelName, Extra, option)
       
        Wb = MCMCPar.Wb;
        Cb = MCMCPar.Cb;
        MeasData = Measurement.MeasData;
        Gamma = MCMCPar.Gamma;

        % initialize ORCHESTRA
%         Extra.ORI = initialize_ORI(Extra.Comp,0);

        evalstr = ['ModPred = ',ModelName,'(x_in,Extra);'];
        eval(evalstr);
        clc
                
        if option == 1, % Model directly computes posterior density
            ptmp = ModPred;
            iitmp = ii;
            
            log_p = log(ptmp);
        end;
        
        if option == 2, % Model computes output simulation
            Err = (MeasData(:)-ModPred(:));
            % Compute the number of measurement data
            N = size(MeasData,1);
            % Derive the log likelihood
            log_p = N.*log(Wb./Sigma) - Cb.*(sum((abs(Err./Sigma)).^(2/(1+Gamma))));
            % And retain in memory
            ptmp = log_p;
                      
        end;
        
        if option == 3, % Model computes output simulation
            Err = (MeasData(:)-ModPred(:));
            % Derive the sum of squared error
            SSR = sum(abs(Err).^(2/(1+Gamma)));
            % And retain in memory
            ptmp = -SSR; log_p = -0.5 * SSR;
        end;
        
        if option == 4, % Model directly computes log posterior density
            ptmp = ModPred;
            log_p = ptmp;
            iitmp = ii;
            log_p = ptmp;
            
        end;
        
        
        if option == 5, % Similar as 3, but now weights with the Measurement Sigma
            % Defime the error
            Err = (MeasData(:)-ModPred(:));
            % Derive the sum of squared error
            SSR = sum(abs(Err).^(2/(1+Gamma)));
            % And retain in memory
            ptmp = -SSR;
            iitmp = ii;
            
            log_p = -0.5 * SSR;
        end;
        
        if option == 6, %SSR weighted with measurement error, maximizing log likelyhood...S.Korteland august 2011
            % Define the error
            Err = (MeasData(:)-ModPred(:));
            % Derive the sum of squared error
            SSR = sum(abs(Err./Sigma).^(2/(1+Gamma)));
            % And retain in memory
            ptmp = -SSR;
            iitmp = ii;
            
            log_p = -0.5 * SSR;
        end;
        
        if option == 8, % Generalized log likelihood (GL)
            % Extract statistical model parameters
            par = Extra.fpar;               % fixed parameters
            par(Extra.idx_vpar) = x_in;  % variable parameters
            par = par';                     % make it a column vector
            statpar = par(end-10:end);
            % Compute the log-likelihood
            log_p = GL('est',statpar,ModPred,MeasData);
            % And retain in memory
            ptmp = log_p;
            iitmp = ii;
            
        end;
    end

end