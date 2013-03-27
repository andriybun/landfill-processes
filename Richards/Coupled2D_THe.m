% Script for running a 2D Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Coupled Simulations for Richard's Equation with Convection Dispersion
% Equation Richard's Equation is solved in a Mixed form with Picards
% Iteration at every spatial difference and Convection Diffusion Equation is
% solved by Marker in Cell method

% Author: Shirishkumar Baviskar and Timo Heimovaara
% Date:7 November 2012

clear
%close all

% Add paths to required functions
addpath('./HelperFunctions');
addpath('./RichardsEquation');
addpath('./DiffusionEquation');
addpath('./MarkerAndCell');

%% Model Geomtery 
% Discretization
% Initialize Solver (discretization of space and time)
% Spatial Discretization
% z and x domain (2d)

s=1;%sides of the domain with origin(0,0) and diagonally opposite corner (s,-s)

% Discretization in x direction
xIN = [0:0.05:0.20];%internodes
nINx = length(xIN);% length of xin internodes
xN(1,1) = xIN(1,1);% leftmost intersection of a node and a internode
xN(1,2:nINx-2) = (xIN(1,2:nINx-2)+xIN(1,3:nINx-1))./2;%middle nodes, avereage of two consecutive internodes is a node
xN(1,nINx-1) = xIN(1,nINx);%right most intersection of a node and a internode
nNx = length(xN);%length of xN nodes
dxN = xN(1,2:end)-xN(1,1:end-1);%differnce in xN nodes
dxIN = xIN(1,2:end)-xIN(1,1:end-1);%differnce in xin inter nodes

% Discretization in z direction
zIN = [-s:0.05:0]';  %internodes in z direction %soil profile until s meters depth
nINz = length(zIN);% length of zIN internodes
zN(1,1) = zIN(1,1);% top intersection of a node and a internode
zN(2:nINz-2,1) = (zIN(2:nINz-2,1)+zIN(3:nINz-1,1))./2;%middle nodes, avereage of two consecutive internodes is a node
zN(nINz-1,1) = zIN(nINz,1);%bottom intersection of a node and a internode
nNz = length(zN);% length of nodes
dzN = zN(2:end,1)-zN(1:end-1,1);% difference in zN nodes
dzIN = zIN(2:end,1)-zIN(1:end-1,1);% difference IN zin internodes



%%
% Model Dimensions
ModelDim.xN=xN;
ModelDim.xIN=xIN;
ModelDim.dxN=dxN;
ModelDim.dxIN=dxIN;
ModelDim.nNx=nNx;
ModelDim.nINx=nINx;

ModelDim.zN=zN;
ModelDim.zIN=zIN;
ModelDim.dzN=dzN;
ModelDim.dzIN=dzIN;
ModelDim.nNz=nNz;
ModelDim.nINz=nINz;

ModelDim.volN = dzIN * dxIN;

%%
% Nodal parameters

allnx = 1:nNx;% 1 to length of all xN nodes
allINx = 1:nINx;% 1 to length of all xin inter nodes
allnz = 1:nNz;% 1 to length of all zN nodes
allINz = 1:nINz;% 1 to length of all zin inter nodes

% Repeated values

ZN=repmat(zN,1,nNx);
XN=repmat(xN,nNz,1);
XXN = repmat(XN,nNz,1);
ZZN = repmat(ZN,1,nNx);
ZZIN = repmat(zIN,1,nINx);
XXIN = repmat(xIN,nINz,1);
ZIN=repmat(zIN,1,nNx);
XIN=repmat(xIN,nNz,1);

%% Richards Equation:
%
% Equation form : Cm.dh/dt=-d(Ksatz.krz)/dz.*(dh/dz+1))-d(Ksatx.krx)/dx.*(dh/dx))
%
% where the dependent variable is h-pressure head
% Cm:Specific mositure capacity
% z: spatial dimension, z negative direction downwards
% x: spatial dimension, x direction parallel to earth surface
% Ksatz: saturated hydraulic conductivity in z direction
% Ksatx: saturated hydraulic conductivity in x direction
% krz: relative permeability in z direction
% krx: relative permeability in x direction
% t: temporal dimesion, time
% Nodal Soil Parameters related to Richards Equation 
SoilPar.alpha(allnz,allnx) = 2.0;% alpha in all nodes in z and x direction
SoilPar.n(allnz,allnx) = 1.5;% n in all all nodes in z and x direction
SoilPar.thetaR(allnz,allnx) = 0.04;% thetaR in all z and x direction
SoilPar.thetaS(allnz,allnx) = 0.4;% thetaS in all z and x direction

% Internodal Soil Parameters
SoilPar.Ksatz (allINz,allnx) = 5e-2; %m/d Saturated hydraulic conductivity in all internodes of z and nodes of x
SoilPar.Ksatx (allnz,allINx) = 5e-2; %m/d Saturated hydraulic conductivity in all nodes of z and intenodes of x

% Impermeable layers analogous to plastics
% Plastics;

%% Boundary Parameters for Richards Equation
% Define boundary condition
% four types of boundary have been implemented: d: Dirichlet, n: Neumann , r: Robbins and z: zero
% flux
% the Neumann is implemented as a flux condition, the Robbins as a leakage condition
%BoundaryPar.BoundaryFunc = @Boundary; % function giving the fluxvalues for all boundary nodes.

% top boundary nodes: Neumann, flux defined in qTop (homogeneous)
iiz = nINz;
iix = allnx;
BoundaryPar.BTypeTB(iiz,iix) = {'n'};
BoundaryPar.NeumannFunc = @qBnd;

% Bottom boundary nodes: Robbins, ambient pressure defined in hAmb (homogeneous)
% external resistance is constant
iiz = 1;
iix = allnx;
BoundaryPar.BTypeTB(iiz,iix) = {'r'};
BoundaryPar.RobbinsFunc = @Ksurf;

% left and right  boundary nodes: Neumann, zero flux condition defined in qTop
iiz = allnz;
iix = [1, nINx];
BoundaryPar.BTypeLR(iiz,iix) = {'n'};


%% Model initial states for Richards Equation;
InitialPar.rhow = 1000; % kg/m3 density of water
InitialPar.g = 9.81; % m/s^2 gravitational acceleration
zref = -1.1; % Height of phreatic surface for initial condition
hbot = -1;% Bottom water head
InitialPar.hIni = hbot+zref-ZN;% initial pressure head in z and x nodes.

% Time for which the Coupled model would be simulated
trange = [0:2.5:50];
% Parameters controlling the Picard Iteration
TimerPar.abstol = 1e-3;
TimerPar.reltol = 1e-3;

TimerPar.dtmax = 10;
% maximum number of iterations before reducing time step
TimerPar.maxiter = 25;
TimerPar.maxiterRed = 25;
TimerPar.succesiter = 10;
% minimum number of iterations for increasing time step
TimerPar.miniter = 15;
% convergence criterium for Picar Iteration, smaller values better
% massbalance however at the cost of more iterations...
% TimerPar.convcrit = 1e-7;
TimerPar.dtini = 1;
TimerPar.IncreaseF = 1.1;
TimerPar.DecreaseF = 0.25;
% please note: time step is also controlled by flux rates...

%% Dynamic Model of Richards Equation

RichardsOutput = SolveRichards(trange,ModelDim,SoilPar,BoundaryPar,TimerPar,InitialPar);
%load('RichardsResults.mat')


%% Convection Dispersion Equation
%
% Equation form :d(theta.C)/dt=-d(theta.Dz)/dz.(dC/dz)+d(theta.vz.C)/dz...
%                            -d(theta.Dx)/dx.(dC/dz)+d(theta.vx.C)/dx
%
% where the dependent variable is C-Concentration
% theta:Moisture of the medium
% z: spatial dimension, z direction below the earth surface
% x: spatial dimension, x direction parallel to earth surface
% Dz: hydrodynamic dispersion in z direction
% Dx: hydrodynamic dispersion in x direction
% vz: adevction velocity in z direction
% vx: adevction velocity in x direction
% t: temporal dimesion, time
%
%% Boundary parameters for Convection Dispersion Equation

% BoundaryPar.CTop = 0;  % Ctop Concentration (Kg/m³),
% 
% %% Model initialization for Convection Dispersion Equation
% % Model Parameters
% 
% Dm=2.03e-9;% Molecular Diffusion of Solute(Cl⁻)in Water [m²/s]
% 
% alphaL=10;%longitudinal dynamic dispersivity [m]
% 
% alphaT=2;%lateral dynamic dispersivity [m]
% %Nodal and Internodal parameters
% 
% %Initial Conditions for CDE
% Ci=1;
% Cini(allnz,allnx)= Ci;
% 
% % Initialization of Markers and Cell approach:
% % number of markers per cell
% nmrkpercell =50;
% ncells = (nINx-1)*(nINz-1);
% nmrk = nmrkpercell*ncells;
% minc = [min(xin(:)),min(zin(:))];
% maxc = [max(xin(:)),max(zin(:))];
% xzmrk = LHSU(minc,maxc,nmrk);
% xmrk = xzmrk(:,1);
% zmrk = xzmrk(:,2);
% nnmrk=length(xzmrk);
% 
% allmrk=1:nnmrk;
% 
% Cmrk(allmrk,1) = Ci;
% 
% %% Dynamic Model of Convection Dispersion Equation
% 
% 
% SoluteTransport2;
% 
% %% PLOTS
% PLOTS;
%% Save Coupled Model
save Coupled2DResults
