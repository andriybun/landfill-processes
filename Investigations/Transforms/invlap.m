% INVLAP  Numerical Inversion of Laplace Transforms
% Source:
% http://www.mathworks.com/matlabcentral/fileexchange/32824-numerical-inversion-of-laplace-transforms-in-matlab/content/INVLAP.m
function [radt, ft] = invlap(FS, tini, tend, nnt, a, ns, nd)
    % Fs is formula for F(s) as a string
    % tini, tend are limits of the solution interval
    % nnt is total number of time instants
    % a, ns, nd are parameters of the method
    % if not given, the method uses implicit values a=6, ns=20, nd=19
    % it is recommended to preserve a=6
    % increasing ns and nd leads to lower error
    % an example of function calling
    % [t,ft]=INVLAP('s/(s^2+4*pi^2)',0,10,1001);
    % to plot the graph of results write plot(t,ft), grid on, zoom on
    tic;
    if nargin == 3
        nnt = 1000;
    end
    if nargin <= 4      % implicit parameters
        a=6;
        ns=20;
        nd=19;
    end
    radt = linspace(tini, tend, nnt); % time vector
    if tini==0          % t=0 is not allowed
        radt = radt(2:1:nnt);
    end
    alpha = zeros(1, ns+1+nd);
    beta = zeros(1, ns+1+nd);
    ft = zeros(1, nnt);
    for n = 1:ns+1+nd   % prepare necessary coefficients
        alpha(n) = a + (n - 1) * pi * 1i;
        beta(n) = -exp(a) * (-1)^n;
    end
    n=1:nd;
    bdif = fliplr(cumsum(gamma(nd+1) ./ gamma(nd+2-n) ./ gamma(n))) ./ 2^nd;
    beta(ns+2:ns+1+nd) = beta(ns+2:ns+1+nd) .* bdif;
    beta(1) = beta(1) / 2;
    for kt = 1:nnt                    % cycle for time t
        tt = radt(kt);
        s = alpha / tt;                % complex frequency s
        bt = beta / tt;
        btF = bt .* FS(s);            % functional value F(s)
        ft(kt) = sum(real(btF));      % original f(tt)
    end
    toc;
end