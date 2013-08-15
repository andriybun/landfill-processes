function [Cb,Wb] = CalcCbWb(Beta);
% This function calculates the parameters for the exponential power density
% Equation [20] paper by Thiemann et al. WRR 2001, Vol 37, No 10, 2521-2535

% First calculate some dummy variables
A1 = gamma(3*(1+Beta)/2); A2 = gamma((1+Beta)/2); 
% And use these to derive Cb and Wb 
Cb = (A1/A2)^(1/(1+Beta)); Wb = sqrt(A1)/((1+Beta)*(A2^(1.5)));