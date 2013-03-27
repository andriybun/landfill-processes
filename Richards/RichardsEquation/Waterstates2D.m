function [States] = Waterstates2D(Time,h,SoilPar,ModelDim,BoundaryPar)
nNx = ModelDim.nNx;
allnx = 1:nNx;
nINx = ModelDim.nINx;
xN=ModelDim.xN;
xIN=ModelDim.xIN;

nNz = ModelDim.nNz;
allnz = 1:nNz;
nINz = ModelDim.nINz;
zN=ModelDim.zN;
zIN=ModelDim.zIN;

h = reshape(h,nNz,nNx);

Zn=repmat(zN,1,nNx);
Xn=repmat(xN,nNz,1);



alpha = SoilPar.alpha;
n = SoilPar.n;
m = 1-1./n;
thetaR = SoilPar.thetaR;
thetaS = SoilPar.thetaS;

Seff = (1+(alpha.*abs(h)).^n).^-m.*(h<0)+(h>=0);
theta = thetaR + Seff.*(thetaS-thetaR);

thetazIN=interp2(Xn,Zn,theta,xN,zIN);% Internodal values of theta i.e
thetaxIN=interp2(Xn,Zn,theta,xIN,zN);% Internodal values of theta i.e
    
[qz,qx] = Richards2D(Time,h,SoilPar,ModelDim,BoundaryPar);
States.Seff = Seff;
States.theta=theta;
States.thetazIN=thetazIN;
States.thetaxIN=thetaxIN;
States.qz=qz;
States.qx=qx;
