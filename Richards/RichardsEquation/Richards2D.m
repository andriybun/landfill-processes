function [qz,qx] = Richards2D(t,h,SoilPar,ModelDim,BoundaryPar)
%Function for calculating rates for a Heat Flow in a spatial domain
%qH are the heat fluxes at all internodal points in the spatial domain;
%t is time
%T are the temperatures in all cells in the spatial domain
%ModelPar contains the Model Parameters (thermal diffusivity of all
%InterNodal points)
%ModelDim contains the ModelDiscretization (coordinates of nodes and 
% internodes, deltas of nodal points and
%internodal points)

nNx = ModelDim.nNx;
allnx = 1:nNx;
nINx = ModelDim.nINx;
xIN = ModelDim.xIN;
xN = ModelDim.xN;

nNz = ModelDim.nNz;
allnz = 1:nNz;
nINz = ModelDim.nINz;
zIN = ModelDim.zIN;
zN = ModelDim.zN;

dxN = repmat(ModelDim.dxN,nNz,1);
dzN = repmat(ModelDim.dzN,1,nNx);

%Calculate Derived States (total head etc.)
[krwinz,krwinx] = RelativePermeabilities(h,SoilPar,ModelDim);
Kx = krwinx.*SoilPar.Ksatx;
Kz = krwinz.*SoilPar.Ksatz;

%% Bottom boundary
% Boundary is either a Neumann (flux) condition or a Robbins (leakage) condition...
Kvalz = zeros(1,nNx);
hAmbz = zeros(1,nNx);
% values for possible Neumann boundary
qvalz = zeros(1,nNx);

% find Neumann conditions
idxn = strcmp(BoundaryPar.BTypeTB(1,allnx),'n');
qtmp = BoundaryPar.NeumannFunc(t,zIN(1),xN(idxn),'bottom',ModelDim);
qvalz(1,idxn) = qtmp;

% find Robbins conditions
idxr = strcmp(BoundaryPar.BTypeTB(1, allnx),'r');
[Ktmp,htmp] = BoundaryPar.RobbinsFunc(t,zIN(1),xN(idxr),'bottom',ModelDim);
Kvalz(1,idxr) = Ktmp;
hAmbz(1,idxr) = htmp;

% please note, vertical nodes are bottom upwards...
qz(1,allnx) = idxn.*qvalz - idxr.*Kvalz.*(h(1,allnx)-hAmbz);

%% middle nodes 
qz(2:nINz-1,allnx)=-Kz(2:nINz-1,allnx).*...
    ((h(2:nNz,allnx)-h(1:nNz-1,allnx))./dzN(1:nNz-1,allnx) + 1);

%% Top boundary
% Boundary is either a Neumann (flux) condition or a Robbins (leakage) condition...
Kvalz = zeros(1,nNx);
hAmbz = zeros(1,nNx);
% values for possible Neumann boundary
qvalz = zeros(1,nNx);

% find Neumann conditions
idxn = strcmp(BoundaryPar.BTypeTB(nINz,allnx),'n');
qtmp = BoundaryPar.NeumannFunc(t,zIN(end),xN(idxn),'top',ModelDim);
qvalz(1,idxn) = qtmp;

% find Robbins conditions
idxr = strcmp(BoundaryPar.BTypeTB(nINz, allnx),'r');
[Ktmp,htmp] = BoundaryPar.RobbinsFunc(t,zIN(end),xN(idxr),'top',ModelDim);
Kvalz(1,idxr) = Ktmp;
hAmbz(1,idxr) = htmp;

qz(nINz,allnx)= idxn.*qvalz - idxr.*Kvalz.*(hAmbz - h(end,allnx));

%% Left boundary
% Boundary is either a Neumann (flux) condition or a Robbins (leakage) condition...
Kvalx = zeros(nNz,1);
hAmbx = zeros(nNz,1);
% values for possible Neumann boundary
qvalx = zeros(nNz,1);

% find Neumann conditions
idxn = strcmp(BoundaryPar.BTypeLR(allnz,1),'n');
qtmp = BoundaryPar.NeumannFunc(t,zN(idxn),xIN(1),'left',ModelDim);
qvalx(idxn,1) = qtmp;

% find Robbins conditions
idxr = strcmp(BoundaryPar.BTypeLR(allnz,1),'r');
[Ktmp,htmp] = BoundaryPar.RobbinsFunc(t,zN(idxr),xIN(1),'left',ModelDim);
Kvalx(idxr,1) = Ktmp;
hAmbx(idxr,1) = htmp;

qx(allnz,1) = idxn.*qtmp - idxr.*Kvalx.*(h(allnz,1)-hAmbx);


%% Middle Nodes
qx(allnz,2:nINx-1)=-Kx(allnz,2:nINx-1).*...
    ((h(allnz,2:nNx)-h(allnz,1:nNx-1))./dxN(allnz,1:nNx-1));% allnz changed from allnx

%% Right boundary
% Boundary is either a Neumann (flux) condition or a Robbins (leakage) condition...
% Boundary is either a Neumann (flux) condition or a Robbins (leakage) condition...
Kvalx = zeros(nNz,1);
hAmbx = zeros(nNz,1);
% values for possible Neumann boundary
qvalx = zeros(nNz,1);

% find Neumann conditions
idxn = strcmp(BoundaryPar.BTypeLR(allnz,nINx),'n');
qtmp = BoundaryPar.NeumannFunc(t,zN(idxn),xIN(end),'right',ModelDim);
qvalx(idxn,1) = qtmp;

% find Robbins conditions
idxr = strcmp(BoundaryPar.BTypeLR(allnz,nINx),'r');
[Ktmp,htmp] = BoundaryPar.RobbinsFunc(t,zN(idxr),xIN(end),'right',ModelDim);
Kvalx(idxr,1) = Ktmp;
hAmbx(idxr,1) = htmp;
qx(allnz,nINx) = idxn.*qtmp - idxr.*Kvalx.*(h(allnz,nNx)-hAmbx);






