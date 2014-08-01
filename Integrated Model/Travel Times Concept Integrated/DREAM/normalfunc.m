function [p] = normalfunc(x,Extra)
% Normal pdf

% Log density
p = -0.5*x*Extra.invC*x';