function S = Kmat(t,h,SoilPar,ModelDim,BoundaryPar)
% Function calculates the jacobian matrix for the Richards equation

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

%Initialize diagonal vectors of Kmat to zero
a = zeros(nNz,nNx);
b = zeros(nNz,nNx);
c = zeros(nNz,nNx);
d = zeros(nNz,nNx);
e = zeros(nNz,nNx);

dxN = repmat(ModelDim.dxN,nNz,1);
dxIN = repmat(ModelDim.dxIN,nINz,1);

dzN = repmat(ModelDim.dzN,1,nNx);
dzIN = repmat(ModelDim.dzIN,1,nINx);

h = reshape(h,nNz,nNx);

%Calculate the relative permeabilities as a function of local water
%pressures
[krwinz,krwinx] = RelativePermeabilities(h,SoilPar,ModelDim);
Kx = krwinx.*SoilPar.Ksatx;
Kz = krwinz.*SoilPar.Ksatz;

%% Middle nodes
iiz = 2:nNz-1;
iix = 2:nNx-1;
a(iiz,iix) = Kx(iiz,iix)./(dxIN(iiz,iix).*dxN(iiz,iix-1));

b(iiz,iix) = Kz(iiz,iix)./(dzIN(iiz,iix).*dzN(iiz-1,iix));

c(iiz,iix) = -(Kz(iiz,iix)./(dzIN(iiz,iix).*dzN(iiz-1,iix)) + ...
   Kz(iiz+1,iix)./(dzIN(iiz,iix).*dzN(iiz,iix)) + ...
   Kx(iiz,iix)./(dxIN(iiz,iix).*dxN(iiz,iix-1)) + ...
   Kx(iiz,iix+1)./(dxIN(iiz,iix).*dxN(iiz,iix)));

d(iiz,iix) = Kz(iiz+1,iix)./(dzIN(iiz,iix).*dzN(iiz,iix));

e(iiz,iix) = Kx(iiz,iix+1)./(dxIN(iiz,iix).*dxN(iiz,iix));%(a,b,c,d,e)

% Boundary Conditions

%Boundary Nodes (corners)
%% Bottom Left Corner
iiz = 1;
iix = 1;
Kvalz = zeros(size(iix));
Kvalx = zeros(size(iix));

%idxn = strcmp(BoundaryPar.BTypeTB(iiz,iix),'n'); 
idxrt = strcmp(BoundaryPar.BTypeTB(iiz,iix),'r');
[Ktmp,~] = BoundaryPar.RobbinsFunc(t,zIN(iiz),xN(idxrt),'bottom',ModelDim);
Kvalz(idxrt) = Ktmp;
idxrl = strcmp(BoundaryPar.BTypeLR(iiz,iix),'r');
[Ktmp,~] = BoundaryPar.RobbinsFunc(t,zN(iiz),xIN(idxrl),'left',ModelDim);
Kvalx(idxrl) = Ktmp;

%a(iiz,iix) = Kx(iiz,iix)./(dxIN(iiz,iix).*dxN(iiz,iix-1));
%b(iiz,iix) = Kz(iiz,iix)./(dzIN(iiz,iix).*dzN(iiz-1,iix));

c(iiz,iix) = -(Kz(iiz+1,iix)./(dzIN(iiz,iix).*dzN(iiz,iix)) + ...
    Kx(iiz,iix+1)./(dxIN(iiz,iix).*dxN(iiz,iix))) - ...
    idxrt.*Kvalz./dzIN(iiz,iix) - idxrl.*Kvalx./dxIN(iiz,iix);

d(iiz,iix) = Kz(iiz+1,iix)./(dzIN(iiz,iix).*dzN(iiz,iix));

e(iiz,iix) = Kx(iiz,iix+1)./(dxIN(iiz,iix).*dxN(iiz,iix));%(c, d, e)


%% Bottom Middle
iiz = 1;
iix = 2:nNx-1;
%idxn = strcmp(BoundaryPar.BTypeTB(iiz,iix),'n'); 
idxr = strcmp(BoundaryPar.BTypeTB(iiz,iix),'r');
Kvalz = zeros(size(idxr));
[Ktmp,~] = BoundaryPar.RobbinsFunc(t,zIN(iiz),xN(idxr),'bottom',ModelDim);
Kvalz(idxr) = Ktmp;

a(iiz,iix) = Kx(iiz,iix)./(dxIN(iiz,iix).*dxN(iiz,iix-1));

%b(iiz,iix) = Kz(iiz,iix)./(dzIN(iiz,iix).*dzN(iiz-1,iix));

c(iiz,iix) = -(Kz(iiz+1,iix)./(dzIN(iiz,iix).*dzN(iiz,iix)) + ...
   Kx(iiz,iix)./(dxIN(iiz,iix).*dxN(iiz,iix-1)) + ...
   Kx(iiz,iix+1)./(dxIN(iiz,iix).*dxN(iiz,iix))) - ...
   idxr.*Kvalz./dzIN(iiz,iix);

d(iiz,iix) = Kz(iiz+1,iix)./(dzIN(iiz,iix).*dzN(iiz,iix));

e(iiz,iix) = Kx(iiz,iix+1)./(dxIN(iiz,iix).*dxN(iiz,iix));%(a,c,d,e)

%% Bottom Right
iiz = 1;
iix = nNx;
Kvalz = zeros(size(iix));
Kvalx = zeros(size(iix));

%idxn = strcmp(BoundaryPar.BTypeTB(iiz,iix),'n'); 
idxrt = strcmp(BoundaryPar.BTypeTB(iiz,iix),'r');
[Ktmp,~] = BoundaryPar.RobbinsFunc(t,zIN(iiz),xN(idxrt),'bottom',ModelDim);
Kvalz(idxrt) = Ktmp;
idxrr = strcmp(BoundaryPar.BTypeLR(iiz,iix+1),'r');
[Ktmp,~] = BoundaryPar.RobbinsFunc(t,zN(iiz),xIN(idxrr),'right',ModelDim);
Kvalx(idxrr) = Ktmp;

a(iiz,iix) = Kx(iiz,iix)./(dxIN(iiz,iix).*dxN(iiz,iix-1));

%b(iiz,iix) = Kz(iiz,iix)./(dzIN(iiz,iix).*dzN(iiz-1,iix));

c(iiz,iix) = -(Kz(iiz+1,iix)./(dzIN(iiz,iix).*dzN(iiz,iix)) + ...
    Kx(iiz,iix)./(dxIN(iiz,iix).*dxN(iiz,iix-1))) - ...
    idxrt.*Kvalz./dzIN(iiz,iix) - idxrr.*Kvalx./dxIN(iiz,iiz);

d(iiz,iix) = Kz(iiz+1,iix)./(dzIN(iiz,iix).*dzN(iiz,iix));

%e(iiz,iix) = Kx(iiz,iix+1)./(dxIN(iiz,iix).*dxN(iiz,iix);%(a,c,d)

%Left Middle 
iiz = 2:nNz-1;
iix = 1;
%idxn = strcmp(BoundaryPar.BTypeLR(iiz,iix),'n'); 
idxr = strcmp(BoundaryPar.BTypeLR(iiz,iix),'r');
Kvalx = zeros(size(idxr));
[Ktmp,~] = BoundaryPar.RobbinsFunc(t,zN(idxr),xIN(iix),'left',ModelDim);
Kvalx(idxr) = Ktmp;
%a(iiz,iix) = Kx(iiz,iix)./(dxIN(iiz,iix).*dxN(iiz,iix-1));

b(iiz,iix) = Kz(iiz,iix)./(dzIN(iiz,iix).*dzN(iiz-1,iix));

c(iiz,iix) = -(Kz(iiz,iix)./(dzIN(iiz,iix).*dzN(iiz-1,iix)) + ...
   Kz(iiz+1,iix)./(dzIN(iiz,iix).*dzN(iiz,iix)) + ...
   Kx(iiz,iix+1)./(dxIN(iiz,iix).*dxN(iiz,iix))) - ...
   idxr.*Kvalx./dxIN(iiz,iix);

d(iiz,iix) = Kz(iiz+1,iix)./(dzIN(iiz,iix).*dzN(iiz,iix));

e(iiz,iix) = Kx(iiz,iix+1)./(dxIN(iiz,iix).*dxN(iiz,iix));%(b,c,d,e)

%Right Middle
iiz = 2:nNz-1;
iix = nNx;
%idxn = strcmp(BoundaryPar.BTypeLR(iiz,iix),'n'); 
idxr = strcmp(BoundaryPar.BTypeLR(iiz,iix+1),'r');
Kvalx = zeros(size(idxr));
[Ktmp,~] = BoundaryPar.RobbinsFunc(t,zN(idxr),xIN(iix+1),'right',ModelDim);
Kvalx(idxr) = Ktmp;

a(iiz,iix) = Kx(iiz,iix)./(dxIN(iiz,iix).*dxN(iiz,iix-1));

b(iiz,iix) = Kz(iiz,iix)./(dzIN(iiz,iix).*dzN(iiz-1,iix));

c(iiz,iix) = -(Kz(iiz,iix)./(dzIN(iiz,iix).*dzN(iiz-1,iix)) + ...
   Kz(iiz+1,iix)./(dzIN(iiz,iix).*dzN(iiz,iix)) + ...
   Kx(iiz,iix)./(dxIN(iiz,iix).*dxN(iiz,iix-1))) - ...
   idxr.*Kvalx./dxIN(iiz,iix);

d(iiz,iix) = Kz(iiz+1,iix)./(dzIN(iiz,iix).*dzN(iiz,iix));

%e(iiz,iix) = Kx(iiz,iix+1)./(dxIN(iiz,iix).*dxN(iiz,iix));%(a,b,c,d)

%% Top Left corner
iiz = nNz;
iix = 1;
Kvalz = zeros(size(iix));
Kvalx = zeros(size(iix));

%idxn = strcmp(BoundaryPar.BTypeTB(iiz,iix),'n'); 
idxrb = strcmp(BoundaryPar.BTypeTB(iiz+1,iix),'r');
[Ktmp,~] = BoundaryPar.RobbinsFunc(t,zIN(iiz+1),xN(idxrb),'top',ModelDim);
Kvalz(idxrb) = Ktmp;
idxrl = strcmp(BoundaryPar.BTypeLR(iiz,iix),'r');
[Ktmp,~] = BoundaryPar.RobbinsFunc(t,zN(iiz),xIN(idxrl),'left',ModelDim);
Kvalx(idxrl) = Ktmp;

%a(iiz,iix) = Kx(iiz,iix)./(dxIN(iiz,iix).*dxN(iiz,iix-1));

b(iiz,iix) = Kz(iiz,iix)./(dzIN(iiz,iix).*dzN(iiz-1,iix));

c(iiz,iix) = -(Kz(iiz,iix)./(dzIN(iiz,iix).*dzN(iiz-1,iix)) + ...
   Kx(iiz,iix+1)./(dxIN(iiz,iix).*dxN(iiz,iix))) - ...
   idxrb.*Kvalz./dzIN(iiz,iix) - idxrl.*Kvalx./dxIN(iiz,iix);

%d(iiz,iix) = Kz(iiz+1,iix)./(dzIN(iiz,iix).*dzN(iiz,iix));

e(iiz,iix) = Kx(iiz,iix+1)./(dxIN(iiz,iix).*dxN(iiz,iix));%(b,c,e)


%Top Middle
iiz = nNz;
iix = 2:nNx-1;
%idxn = strcmp(BoundaryPar.BTypeTB(iiz+1,iix),'n'); 
idxr = strcmp(BoundaryPar.BTypeTB(iiz+1,iix),'r');
Kvalz = zeros(size(idxr));
[Ktmp,~] = BoundaryPar.RobbinsFunc(t,zIN(iiz+1),xN(idxr),'top',ModelDim);
Kvalz(idxr) = Ktmp;

a(iiz,iix) = Kx(iiz,iix)./(dxIN(iiz,iix).*dxN(iiz,iix-1));

b(iiz,iix) = Kz(iiz,iix)./(dzIN(iiz,iix).*dzN(iiz-1,iix));

c(iiz,iix) = -(Kz(iiz,iix)./(dzIN(iiz,iix).*dzN(iiz-1,iix)) + ...
   Kx(iiz,iix)./(dxIN(iiz,iix).*dxN(iiz,iix-1)) + ...
   Kx(iiz,iix+1)./(dxIN(iiz,iix).*dxN(iiz,iix))) - ...
   idxrb.*Kvalz./dzIN(iiz,iix);

%d(iiz,iix) = Kz(iiz+1,iix)./(dzIN(iiz,iix).*dzN(iiz,iix));

e(iiz,iix) = Kx(iiz,iix+1)./(dxIN(iiz,iix).*dxN(iiz,iix));%(a,b,c,e)

%Top Right
iiz = nNz;
iix = nNx;
Kvalz = zeros(size(iix));
Kvalx = zeros(size(iix));

%idxn = strcmp(BoundaryPar.BTypeTB(iiz,iix),'n'); 
idxrb = strcmp(BoundaryPar.BTypeTB(iiz+1,iix),'r');
[Ktmp,~] = BoundaryPar.RobbinsFunc(t,zIN(iiz+1),xN(idxrb),'top',ModelDim);
Kvalz(idxrb) = Ktmp;
idxrr = strcmp(BoundaryPar.BTypeLR(iiz,iix+1),'r');
[Ktmp,~] = BoundaryPar.RobbinsFunc(t,zN(iiz),xIN(idxrr),'right',ModelDim);
Kvalx(idxrr) = Ktmp;

a(iiz,iix) = Kx(iiz,iix)./(dxIN(iiz,iix).*dxN(iiz,iix-1));

b(iiz,iix) = Kz(iiz,iix)./(dzIN(iiz,iix).*dzN(iiz-1,iix));

c(iiz,iix) = -(Kz(iiz,iix)./(dzIN(iiz,iix).*dzN(iiz-1,iix)) + ...
   Kx(iiz,iix)./(dxIN(iiz,iix).*dxN(iiz,iix-1))) - ...
   idxrb.*Kvalz./dzIN(iiz,iix) - idxrr.*Kvalx./dxIN(iiz,iix);

%d(iiz,iix) = Kz(iiz+1,iix)./(dzIN(iiz,iix).*dzN(iiz,iix));

%e(iiz,iix) = Kx(iiz,iix+1)./(dxIN(iiz,iix).*dxN(iiz,iix);%(a,b,c);



B =  diag(a(nNz+1:end),-nNz) + diag(b(2:end),-1)+diag(c(:),0)+...
    diag(d(1:end-1),+1) + diag(e(1:end-nNz),nNz);
S = sparse(B);
end