function [krwinz,krwinx] =RelativePermeabilities(h,SoilPar,ModelDim)
nNx = ModelDim.nNx;% length of xn
nNz = ModelDim.nNz;% length of zn
allnx = 1:nNx;% 1 to xn
nINx = ModelDim.nINx;% length of xin 
allnz = 1:nNz;%1 to zn
nINz = ModelDim.nINz;% length of zin

h = reshape(h,nNz,nNx);
alpha = SoilPar.alpha(allnz,allnx);% alpha in z and x nodes
n = SoilPar.n(allnz,allnx);% n in x and z nodes
m = 1-1./n;% m in z an x nodes

Seff = (1+(alpha.*abs(h)).^n).^-m.*(h<0)+(h>=0);% Effective Saturation
krw = Seff.^0.5 .* (1-(1-Seff.^(1./m)).^m).^2.*(h<0)+(h>=0);% relative permeability

krwtmp = reshape(krw,nNz,nNx);% relative permeability in z and x nodal direction
%% to do Estimate kr for internodes...
krwinz(1,allnx) = krwtmp(1,allnx);
krwinz(2:nINz-1,allnx) = min (krwtmp(1:nNz-1,allnx),krwtmp(2:nNz,allnx));
krwinz(nINz,allnx) = krwtmp(nNz,allnx);

krwinx(allnz,1) = krwtmp(allnz,1);
krwinx(allnz,2:nINx-1) = min (krwtmp(allnz,1:nNx-1),krwtmp(allnz,2:nNx));
krwinx(allnz,nINx) = krwtmp(allnz,nNx);
