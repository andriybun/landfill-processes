function [newgen,accept] = metrop(x,p_x,log_p_x,x_old,p_old,log_p_old,alpha_s,Measurement,MCMCPar,Extra,option);
% Metropolis rule for acceptance or rejection

% Calculate the number of Chains
NrChains = size(x,1);

% First set newgen to the old positions in X
newgen = [x_old p_old log_p_old];

% And initialize accept with zeros
accept = zeros(NrChains,1);

% -------------------- Now check the various options ----------------------
if option == 1,
    alpha = p_x./p_old;
end;

if option == 2 | option == 4 | option == 8, % Lnp probability evaluation
    alpha = exp(p_x - p_old);
end;

if option == 3, % SSE probability evaluation
    alpha = (p_x./p_old).^(-Measurement.N.*(1+MCMCPar.Gamma)./2);
end;

if option == 5, % Similar as 3 but now weighted with Measurement.Sigma
    alpha = exp(-0.5*(-p_x + p_old)./Measurement.Sigma^2); % signs are different because we write -SSR
end;

if option == 6, %S.Korteland august 2011
    alpha = exp(0.5*(p_x - p_old)); % 
end;
% -------------------------------------------------------------------------

% Using a prior or not? If so, alpha needs to adjusted: multiplied with prior
% -------------------------------------------------------------------------
if strcmp(Extra.InitPopulation,'PRIOR');
	
    % Compute prior densities for each parameter in each sequence
    for qq = 1:MCMCPar.n,
        for zz = 1:MCMCPar.seq,
            % Compute prior of proposal 
            prior_x(zz,qq) = eval(char(strrep(Extra.prior(qq),'rnd(','pdf(x(zz,qq),')));
            % Compute prior of current location
            prior_old(zz,qq) = eval(char(strrep(Extra.prior(qq),'rnd(','pdf(x_old(zz,qq),')));
        end;
    end;

    % Take the product of the various densities
    prior_old = max(prod(prior_old,2),1e-299); prior_prop = max(prod(prior_x,2),1e-299); % (to avoid 0/0 --> NaN)
    % Take the ratio
    alpha_pr = prior_prop./prior_old;
    % Now update alpha value with prior
    alpha = alpha.*alpha_pr;  

end;
% -------------------------------------------------------------------------

% Modify 
alpha = alpha.* alpha_s;

% Generate random numbers
Z = rand(NrChains,1);

% Find which alpha's are greater than Z
idx = find(alpha > Z);

% And update these chains
newgen(idx,1:MCMCPar.n+2) = [x(idx,1:MCMCPar.n) p_x(idx,1) log_p_x(idx,1)];

% And indicate that these chains have been accepted
accept(idx,1) = 1;