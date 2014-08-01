function [p] = mixturemodel(x,Extra);
% Computes SSR of bananafunc

% Forward call
p = 1/3 * mvnpdf(x,Extra.mu1,Extra.sigma) + 2/3 * mvnpdf(x,Extra.mu2,Extra.sigma);

p = diag(p);
