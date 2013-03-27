function S = Yvec(t,h,SoilPar,ModelDim,BoundaryPar)
% Function calculates the RHS the Matrix Solution of the Richards equation

nNz = ModelDim.nNz;
allnz = 1:nNz;
nINz = ModelDim.nINz;
allinz = 1:nINz;
zIN = ModelDim.zIN;
zN = ModelDim.zN;

nNx = ModelDim.nNx;
allnx = 1:nNx;
nINx = ModelDim.nINx;
allinx = 1:nINx;
xIN = ModelDim.xIN;
xN = ModelDim.xN;

dxN = repmat(ModelDim.dxN,nNz,1);
dxIN = repmat(ModelDim.dxIN,nINz,1);

dzN = repmat(ModelDim.dzN,1,nNx);
dzIN = repmat(ModelDim.dzIN,1,nINx);

h = reshape(h,nNz,nNx);

% Calculate the relative permeabilities as a function of local water
% pressures
[krwinz,krwinx] = RelativePermeabilities(h,SoilPar,ModelDim);
Kx = krwinx.*SoilPar.Ksatx;
Kz = krwinz.*SoilPar.Ksatz;


%% Middle nodes
% Only gravity terms... 
iiz = 2:nNz-1;
iix = 2:nNx-1;

Y(iiz,iix) = (Kz(iiz+1,iix)-Kz(iiz,iix))./dzIN(iiz,iix);

% Boundary Nodes:
% Two types of boundary conditions are implemented: Neumann and Robbins
% The boundary values are given by the functions passed with BoudaryPar.
% the Neumann boundary function passes the boundary flux, the Robbins boundary condition passes the
% resistance value and the ambient boundary pressure...

% Boundary Nodes (corners)
%% Bottom left
iiz = 1;
iix = 1;
% values for possible Robbins boundary
Kvalz = zeros(length(iix));
Kvalx = zeros(length(iix));
hAmbz = zeros(length(iix));
hAmbx = zeros(length(iix));
% values for possible Neumann boundary
qvalx = zeros(length(iix));
qvalz = zeros(length(iix));

idxnt = strcmp(BoundaryPar.BTypeTB(iiz,iix),'n');
qtmp = BoundaryPar.NeumannFunc(t,zIN(iiz),xN(idxnt),'bottom',ModelDim);
qvalz(idxnt) = qtmp;

idxnl = strcmp(BoundaryPar.BTypeLR(iiz,iix),'n');
qtmp = BoundaryPar.NeumannFunc(t,zN(iiz),xIN(idxnl),'left',ModelDim);
qvalx(idxnl) = qtmp;

idxrt = strcmp(BoundaryPar.BTypeTB(iiz,iix),'r');
[Ktmp,htmp] = BoundaryPar.RobbinsFunc(t,zIN(iiz),xN(idxrt),'bottom',ModelDim);
Kvalz(idxrt) = Ktmp;
hAmbz(idxrt) = htmp;

idxrl = strcmp(BoundaryPar.BTypeLR(iiz,iix),'r');
[Ktmp,htmp] = BoundaryPar.RobbinsFunc(t,zN(iiz),xIN(idxrl),'left',ModelDim);
Kvalx(idxrl) = Ktmp;
hAmbx(idxrl) = htmp;

Y(iiz,iix) = Kz(iiz+1,iix)./dzIN(iiz,iix) + ...
    idxnt.*qvalz./dzIN(iiz,iix) + idxnl.*qvalx./dxIN(iiz,iix) + ...
    idxrt.*Kvalz.*hAmbz./dzIN(iiz,iix) + idxrl.*Kvalx.*hAmbx./dxIN(iiz,iix);

%% Bottom middle
iiz = 1;
iix = 2:nNx-1;
% values for possible Robbins boundary
Kvalz = zeros(1,length(iix));
hAmbz = zeros(1,length(iix));
% values for possible Neumann boundary
qvalz = zeros(1,length(iix));

idxnt = strcmp(BoundaryPar.BTypeTB(iiz,iix),'n');
qtmp = BoundaryPar.NeumannFunc(t,zIN(iiz),xN(idxnt),'bottom',ModelDim);
qvalz(1,idxnt) = qtmp;

idxrt = strcmp(BoundaryPar.BTypeTB(iiz,iix),'r');
[Ktmp,htmp] = BoundaryPar.RobbinsFunc(t,zIN(iiz),xN(idxrt),'bottom',ModelDim);
Kvalz(1,idxrt) = Ktmp;
hAmbz(1,idxrt) = htmp;


Y(iiz,iix) = Kz(iiz+1,iix)./dzIN(iiz,iix) + ...
    idxnt.*qvalz./dzIN(iiz,iix) + ...
    idxrt.*Kvalz.*hAmbz./dzIN(iiz,iix);

%% Bottom Right
iiz = 1;
iix = nNx;
% values for possible Robbins boundary
Kvalz = zeros(length(iix));
Kvalx = zeros(length(iix));
hAmbz = zeros(length(iix));
hAmbx = zeros(length(iix));
% values for possible Neumann boundary
qvalx = zeros(length(iix));
qvalz = zeros(length(iix));

idxnt = strcmp(BoundaryPar.BTypeTB(iiz,iix),'n');
qtmp = BoundaryPar.NeumannFunc(t,zIN(iiz),xN(idxnt),'bottom',ModelDim);
qvalz(idxnt) = qtmp;

idxnr = strcmp(BoundaryPar.BTypeLR(iiz,iix+1),'n');
qtmp = BoundaryPar.NeumannFunc(t,zN(iiz),xIN(idxnr),'right',ModelDim);
qvalx(idxnr) = qtmp;

idxrt = strcmp(BoundaryPar.BTypeTB(iiz,iix),'r');
[Ktmp,htmp] = BoundaryPar.RobbinsFunc(t,zIN(iiz),xN(idxrt),'bottom',ModelDim);
Kvalz(idxrt) = Ktmp;
hAmbz(idxrt) = htmp;

idxrr = strcmp(BoundaryPar.BTypeLR(iiz,iix+1),'r');
[Ktmp,htmp] = BoundaryPar.RobbinsFunc(t,zN(iiz),xIN(idxrr),'right',ModelDim);
Kvalx(idxrr) = Ktmp;
hAmbx(idxrr) = htmp;

Y(iiz,iix) = Kz(iiz+1,iix)./dzIN(iiz,iix) + ...
    idxnt.*qvalz./dzIN(iiz,iix) - idxnr.*qvalx./dxIN(iiz,iix) + ...
    idxrt.*Kvalz.*hAmbz./dzIN(iiz,iix) + idxrr.*Kvalx.*hAmbx./dxIN(iiz,iix);

%% Left Middle 
iiz = 2:nNz-1;
iix = 1;
% values for possible Robbins boundary
Kvalx = zeros(length(iiz),1);
hAmbx = zeros(length(iiz),1);
% values for possible Neumann boundary
qvalx = zeros(length(iiz),1);

idxnl = strcmp(BoundaryPar.BTypeLR(iiz,iix),'n');
qtmp = BoundaryPar.NeumannFunc(t,zN(idxnl),xIN(iix),'left',ModelDim);
qvalx(idxnl,1) = qtmp;

idxrl = strcmp(BoundaryPar.BTypeLR(iiz,iix),'r');
[Ktmp,htmp] = BoundaryPar.RobbinsFunc(t,zN(idxrl),xIN(iix),'left',ModelDim);
Kvalx(idxrl,1) = Ktmp;
hAmbx(idxrl,1) = htmp;

Y(iiz,iix) = Kz(iiz+1,iix)./dzIN(iiz,iix)-Kz(iiz,iix)./dzIN(iiz,iix) + ...
    idxnl.*qvalx./dxIN(iiz,iix) + ...
    idxrl.*Kvalx.*hAmbx./dxIN(iiz,iix);

%% Right Middle
iiz = 2:nNz-1;
iix = nNx;
% values for possible Robbins boundary
Kvalx = zeros(length(iiz),1);
hAmbx = zeros(length(iiz),1);
% values for possible Neumann boundary
qvalx = zeros(length(iiz),1);

idxnr = strcmp(BoundaryPar.BTypeLR(iiz,iix+1),'n');
qtmp = BoundaryPar.NeumannFunc(t,zN(idxnr),xIN(iix+1),'right',ModelDim);
qvalx(idxnr,1) = qtmp;

idxrr = strcmp(BoundaryPar.BTypeLR(iiz,iix+1),'r');
[Ktmp,htmp] = BoundaryPar.RobbinsFunc(t,zN(idxrr),xIN(iix+1),'right',ModelDim);
Kvalx(idxrr,1) = Ktmp;
hAmbx(idxrr,1) = htmp;

Y(iiz,iix) = Kz(iiz+1,iix)./dzIN(iiz,iix)-Kz(iiz,iix)./dzIN(iiz,iix) - ...
    idxnr.*qvalx./dxIN(iiz,iix) + ...
    idxrr.*Kvalx.*hAmbx./dxIN(iiz,iix);

%% Top Left Corner
iiz = nNz;
iix = 1;
Kvalz = zeros(length(iix));
Kvalx = zeros(length(iix));
hAmbz = zeros(length(iix));
hAmbx = zeros(length(iix));
% values for possible Neumann boundary
qvalx = zeros(length(iix));
qvalz = zeros(length(iix));

idxnb = strcmp(BoundaryPar.BTypeTB(iiz+1,iix),'n');
qtmp = BoundaryPar.NeumannFunc(t,zIN(iiz+1),xN(idxnb),'top',ModelDim);
qvalz(idxnb) = qtmp;

idxnl = strcmp(BoundaryPar.BTypeLR(iiz,iix),'n');
qtmp = BoundaryPar.NeumannFunc(t,zN(iiz),xIN(idxnl),'left',ModelDim);
qvalx(idxnl) = qtmp;

idxrb = strcmp(BoundaryPar.BTypeTB(iiz+1,iix),'r');
[Ktmp,htmp] = BoundaryPar.RobbinsFunc(t,zIN(iiz+1),xN(idxrb),'top',ModelDim);
Kvalz(idxrb) = Ktmp;
hAmbz(idxrb) = htmp;

idxrl = strcmp(BoundaryPar.BTypeLR(iiz,iix),'r');
[Ktmp,htmp] = BoundaryPar.RobbinsFunc(t,zN(iiz),xIN(idxrl),'left',ModelDim);
Kvalx(idxrl) = Ktmp;
hAmbx(idxrl) = htmp;

Y(iiz,iix) = -Kz(iiz,iix)./dzIN(iiz,iix) - ...
    idxnb.*qvalz./dzIN(iiz,iix) + idxnl.*qvalx./dxIN(iiz,iix) + ...
    idxrb.*Kvalz.*hAmbz./dzIN(iiz,iix) + idxrl.*Kvalx.*hAmbx./dxIN(iiz,iix);

%% Top Middle
iiz = nNz;
iix = 2:nNx-1;
Kvalz = zeros(1,length(iix));
hAmbz = zeros(1,length(iix));
% values for possible Neumann boundary
qvalz = zeros(1,length(iix));

idxnb = strcmp(BoundaryPar.BTypeTB(iiz+1,iix),'n');
qtmp = BoundaryPar.NeumannFunc(t,zIN(iiz+1),xN(idxnb),'top',ModelDim);
qvalz(1,idxnb) = qtmp;

idxrb = strcmp(BoundaryPar.BTypeTB(iiz+1,iix),'r');
[Ktmp,htmp] = BoundaryPar.RobbinsFunc(t,zIN(iiz+1),xN(idxrb),'top',ModelDim);
Kvalz(1,idxrb) = Ktmp;
hAmbz(1,idxrb) = htmp;

Y(iiz,iix) = -Kz(iiz,iix)./dzIN(iiz,iix) - ...
    idxnb.*qvalz./dzIN(iiz,iix) + ...
    idxrb.*Kvalz.*hAmbz./dzIN(iiz,iix);

%% Top Right
iiz = nNz;
iix = nNx;
Kvalz = zeros(length(iix));
Kvalx = zeros(length(iix));
hAmbz = zeros(length(iix));
hAmbx = zeros(length(iix));
% values for possible Neumann boundary
qvalx = zeros(length(iix));
qvalz = zeros(length(iix));

idxnb = strcmp(BoundaryPar.BTypeTB(iiz+1,iix),'n');
qtmp = BoundaryPar.NeumannFunc(t,zIN(iiz+1),xN(idxnb),'top',ModelDim);
qvalz(idxnb) = qtmp;

idxnr = strcmp(BoundaryPar.BTypeLR(iiz,iix+1),'n');
qtmp = BoundaryPar.NeumannFunc(t,zN(iiz),xIN(iix+1),'right',ModelDim);
qvalx(idxnr) = qtmp;

idxrb = strcmp(BoundaryPar.BTypeTB(iiz+1,iix),'r');
[Ktmp,htmp] = BoundaryPar.RobbinsFunc(t,zIN(iiz+1),xN(iix),'top',ModelDim);
Kvalz(idxrb) = Ktmp;
hAmbz(idxrb) = htmp;

idxrr = strcmp(BoundaryPar.BTypeLR(iiz,iix+1),'r');
[Ktmp,htmp] = BoundaryPar.RobbinsFunc(t,zN(iiz),xIN(idxrr),'right',ModelDim);
Kvalx(idxrr) = Ktmp;
hAmbx(idxrr) = htmp;

Y(iiz,iix) = -Kz(iiz,iix)./dzIN(iiz,iix) - ...
    idxnb.*qvalz./dzIN(iiz,iix) - idxnr.*qvalx./dxIN(iiz,iix) + ...
    idxrb.*Kvalz.*hAmbz./dzIN(iiz,iix) + idxrr.*Kvalx.*hAmbx./dxIN(iiz,iix);

B =  Y;
S = B(:);
end 