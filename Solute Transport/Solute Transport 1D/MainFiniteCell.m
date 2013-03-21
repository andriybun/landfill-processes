%% Main script for running Finite Cell Method (Sun N.-Z., 1999) in combination with Marker-in-Cell
%  to simulate advection-dispersion processes

clc
clear

% Path to directory with code for solving Richards equation
addpath('Common');
addpath('PhysicalProcesses');
addpath('PhysicalProcesses/Richards');

% Numerical tolerance
SimulationPar = struct();
SimulationPar.EPSILON = 1e-10;

%% Model parameters
%Initialize Solver (discretization of space and time)

% Spatial Discretization
% Define internodal boundaries (include at least all soillayer boundaries)
zTop = 0;
zBottom = -1;
dz = -0.05;
ModelDim = InitializeNodes('z', zTop, zBottom, dz);

nSolutes = 1;

alln = 1:ModelDim.znn;
allin = 1:ModelDim.znin;

% Boundary parameters
BoundaryPar.cTop = zeros(1, nSolutes);
BoundaryPar.qTop = @qBoundary;
BoundaryPar.kSurf = 1e-2;
BoundaryPar.hAmb = -0.025;
BoundaryPar.zRef = 0.2;        % Height of phreatic surface for initial condition
BoundaryPar.hBot = zBottom;

% Time discretization
tEnd = 200;
dt = 2;
tRange = 0:dt:tEnd;
t = tRange(1);

SimulationPar.dt = dt;
SimulationPar.tRange = tRange;
SimulationPar.TIME_EPSILON = 1e-7;

kSat = 1e-2;                             % m^3 water per m^2 soil per day

% Hydraulic properties
SoilPar = InitializeSoilProperties(kSat, ModelDim);

% Diffusion coefficients
SoilPar.d = 0 * 1e-3; %[1e-3, 1e-4];

%% 
tic

% Max timestep is based on cell size and flow velocity
dtMax = max(abs(ModelDim.dzn)) ./ max(abs(kSat));

% Initialize output matrices
timeOutVec = nan(size(tRange));
timeOutVec(1) = t;
iTime = 1;

%% Solve Richards
tic
[qOutR, thetaOutR, hOutR, tRangeR] = ...
    Richards(tRange, dtMax, ModelDim, SoilPar, BoundaryPar);


%% Validation
validationMode = false;
if validationMode
    qOutR = -1e-2 * ones(size(qOutR));
    thetaOutR = SoilPar.thetaS * ones(size(thetaOutR));
    inFlow = min(min(qOutR)) / SoilPar.thetaS;
    doDisplayAnalyticalSolution = true;
else
    doDisplayAnalyticalSolution = false;
    inFlow = min(BoundaryPar.qTop(tRange)) / SoilPar.thetaS;
end
toc

%% Initialize markers object
MarkerData = MarkerDataCl(thetaOutR(1, :)', nSolutes, ...
    ModelDim, SoilPar, SimulationPar, @InitialConcentration);

%% Initialize outputs
initialMass = MarkerData.dv' * MarkerData.c;
mSoluteOutCum = zeros(1, nSolutes);
cOutArr = nan(ModelDim.znn, nSolutes, numel(tRange));
qOut = zeros(1, numel(tRange));
qOut(1) = -qOutR(1, ModelDim.znin);
mOut = zeros(nSolutes, numel(tRange));

[MarkerData, cOutArr(:, :, iTime)] = MarkerData.NodalConcentrations();
mOut(:, 1) = qOut(1) * cOutArr(ModelDim.znn, :, iTime);

%% Simulation loop
tic

thetaNextN = interp2(ModelDim.zn, tRangeR, thetaOutR, ModelDim.zn, t)';
isMidStep = false;

while abs(t - tEnd) > SimulationPar.TIME_EPSILON,
    dtRemaining = min(tRange(iTime + 1) - t, SimulationPar.dt);

    % Recalculate all that stuff only if proceeding to a next time step
    if (~isMidStep)
        tPrev = t;
        vOutStep = 0;
        mSoluteOutStep = 0;

        % Interpolate values of moisture content for both ends of time Interval
        thetaN = thetaNextN;
        thetaNextN = interp2(ModelDim.zn, tRangeR, thetaOutR, ModelDim.zn, t + dtRemaining)';
        
        % Calculate mass preserving internodal flux for given time step
        deltaQn = ModelDim.dzin ./ dtRemaining .* (thetaNextN - thetaN);
        qIn = zeros(ModelDim.znin, 1);
        qIn(1) = interp1(tRangeR, qOutR(:, 1), t + dtRemaining);
        qIn(2:ModelDim.znin) = qIn(1) - cumsum(deltaQn);
    end
    
    % Advect markers, inject new markers, recompute concentrations
    [MarkerData, dtCalc, vOut, mSoluteOut] = ...
        MarkerData.Advect(t, dtRemaining, qIn, thetaN, BoundaryPar);
    [MarkerData, thetaNewN] = MarkerData.NodalThetas();

    % Accumulate discharge for current time step
    vOutStep = vOutStep + vOut;
    
    % Mass balance and transfer correctness check
    mSoluteOutStep = mSoluteOutStep + mSoluteOut;
    mSoluteOutCum = mSoluteOutCum + mSoluteOut;
    
    % Inclrease time
    t = t + dtCalc;
    
    % Mass balance and transfer correctness check: expected theta
    thetaExp = thetaN - (t - tPrev) * (qIn(2:end) - qIn(1:end-1)) ./ ModelDim.dzin;
    checkExp = thetaNewN - thetaExp;
    if (any(~RealEq(checkExp, 0, SimulationPar.EPSILON)))
        warning('ResultValidation:MassBalance', 't = %5.3f: Fluid mass balance error: %5.3e', ...
            t - dtCalc, max(abs(checkExp)));
    end
    mRemaining = MarkerData.dv' * MarkerData.c;
    balance = abs(initialMass - mSoluteOutCum - mRemaining);
    if (abs(balance) > SimulationPar.EPSILON)
        warning('ResultValidation:MassBalance', 't = %5.3f: Solute mass balance error: %5.3e', ...
            t - dtCalc, abs(balance));
    end

    % Check if time step was reduced during advection
    if (~RealEq(dtCalc, dtRemaining, SimulationPar.EPSILON))
        isMidStep = true;
    else
        % Update output matrices
        iTime = iTime + 1;
        [MarkerData, cOutArr(:, :, iTime)] = MarkerData.NodalConcentrations();
        timeOutVec(iTime) = t;
        qOut(iTime) = vOutStep;
        mOut(:, iTime) = mSoluteOutStep;
        
        % Reset flag
        isMidStep = false;
    end
end
toc

tic
%% Running analytical solution
cAnalyticalArr = SoluteTransportAnalytic(ModelDim.zn, timeOutVec, inFlow, SoilPar.d(1), ModelDim);
toc

%% Checking mass balance
fprintf('\nMass balance check:\n');
mRemaining = MarkerData.dv' * MarkerData.c;
balance = abs(initialMass - mSoluteOutCum - mRemaining);
for soluteIdx = 1:nSolutes
    fprintf('Solute %d: ', soluteIdx);
    fprintf('\tInitial mass: %f, flushed %f, remaining %f\n', ...
        initialMass(soluteIdx), mSoluteOutCum(soluteIdx), mRemaining(soluteIdx));
    fprintf('\t\tBalance is: %e\n', balance(soluteIdx));
    if (abs(balance(soluteIdx)) > SimulationPar.EPSILON)
        warning('ResultValidation:MassBalance', 'Final check: solute mass balance error!');
    end
end
    
%% Plots
if zTop == 3
    selectedNodes = [20:20:60, 70, 80, 100:20:ModelDim.znn];
elseif zTop == 0
    if zBottom == -1
        selectedNodes = [5:5:20];
    else
        if SoilPar.d == 0
            selectedNodes = [1, 15, 25, 50:25:ModelDim.znn];
        else
            selectedNodes = [1, 10, 20, 40:20:ModelDim.znn];
        end
    end
else
    firstNodeIdx = 5;
    displayStep = 5;
    selectedNodes = firstNodeIdx:displayStep:ModelDim.znn;
end
DisplayPlotSinglePhase(timeOutVec, cOutArr, cAnalyticalArr, ModelDim, SoilPar, inFlow, ...
    selectedNodes, doDisplayAnalyticalSolution);

% Plot outflux vs. out concentration
figure(2)
clf
plotyy(timeOutVec, qOut, timeOutVec, mOut(soluteIdx, :) ./ qOut);