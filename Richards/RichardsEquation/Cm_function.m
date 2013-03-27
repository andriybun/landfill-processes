function S = Cm_function(t,h,SoilPar,ModelDim)
nNx = ModelDim.nNx;
allnx = 1:nNx;
nINx = ModelDim.nINx;

nNz = ModelDim.nNz;
allnz = 1:nNz;
nINz = ModelDim.nINz;

alpha = SoilPar.alpha(allnz,allnx);
n = SoilPar.n(allnz,allnx);
m = 1-1./n;

Seff = (1+(alpha.*abs(h)).^n).^-m.*(h<0)+(h>=0);
S = Seff+SoilPar.thetaR./SoilPar.thetaS;

C = m.*n.*alpha.*(alpha.*abs(h)).^(n-1).*(1+(alpha.*abs(h)).^n).^(-m-1) .* ...
   (h<0)+0.*(h>=0);
Sw = 4e-10.*1000.*9.81; % compressibility of water [1/m]
C = C+Sw.*S;
S = sparse(diag(C(:),0));