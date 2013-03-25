%% Main script for running Finite Cell Method (Sun N.-Z., 1999) in combination with Marker-in-Cell
%  to simulate advection-dispersion processes integrated with diffusion processes between
%  mobile and immobile phases.

clc
clear

% Path to directory with code for solving Richards equation
addpath('../../Common');
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

diffCoeff = 1e-3; % [1e-3, 5e-3];
exchCoeff = 0 * 3e-3; % [1e-1, 5e-1];
nSolutes = numel(diffCoeff);

alln = 1:ModelDim.znn;
allin = 1:ModelDim.znin;

% Boundary parameters
BoundaryPar.cTop = zeros(1, nSolutes);
BoundaryPar.qTop = @qBoundary;
BoundaryPar.kSurf = 1e-2;
BoundaryPar.hAmb = -0.025;
BoundaryPar.zRef = 0;       % Height of phreatic surface for initial condition
BoundaryPar.hBot = zBottom;

% Time discretization
tEnd = 300;
dt = 2;
tRange = 0:dt:tEnd;
t = tRange(1);

SimulationPar.dt = dt;
SimulationPar.tRange = tRange;
SimulationPar.TIME_EPSILON = 1e-7;

% Initialize fraction of mobile phase volume
% ModelDim.mobileFraction = linspace(0.9, 0.1, ModelDim.znn)';
ModelDim.mobileFraction = 0.1 * ones(ModelDim.znn, 1);

% Saturated conductivity
kSat = 1e-2;             % m^3 water per m^2 soil per day

% Hydraulic properties
SoilPar = InitializeSoilProperties(kSat, ModelDim); % , ModelDim.mobileFraction

% Diffusion coefficients
SoilPar.d = diffCoeff;
% Exchange rates between phases
SoilPar.kExch = exchCoeff;

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
inFlow = min(BoundaryPar.qTop(tRange));
toc

%% Initialize markers object
MarkerData = MarkerDataCl(thetaOutR(1, :)', nSolutes, ...
    ModelDim, SoilPar, SimulationPar, @InitialConcentration);

%% Initialize immoble phase object
thetaN = interp2(ModelDim.zn, tRangeR, thetaOutR, ModelDim.zn, t)';
ImmobilePhase = ImmobilePhaseCl(thetaN, nSolutes, ModelDim, SoilPar, SimulationPar);

%% Initialize some variables for checking mass balance
mInitialMobile = MarkerData.dv' * MarkerData.c;
mInitialImmobile = (-ModelDim.dzin .* ImmobilePhase.theta)' * ImmobilePhase.c;
mInitial = mInitialMobile + mInitialImmobile;

%% Initialize outputs
mSoluteOutCum = zeros(1, nSolutes);
cMobOutArr = nan(ModelDim.znn, nSolutes, numel(tRange));
cImmobOutArr = nan(ModelDim.znn, nSolutes, numel(tRange));
qOut = zeros(1, numel(tRange));
qOut(1) = -qOutR(1, ModelDim.znin);
mOut = zeros(nSolutes, numel(tRange));

[MarkerData, cMobOutArr(:, :, iTime)] = MarkerData.NodalConcentrations();
mOut(:, 1) = qOut(1) * cMobOutArr(ModelDim.znn, :, iTime);
cImmobOutArr(:, :, iTime) = ImmobilePhase.c;

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
    
    % Compute concentrations of solutes in mobile phase
    [MarkerData, cMob] = MarkerData.NodalConcentrations();
    
    % Compute diffusion between mobile and mobile phase
    [ImmobilePhase, cMob] = ImmobilePhase.Diffuse(t, dtRemaining, cMob, thetaNewN);
    MarkerData = MarkerData.ApplySubgridDiffusion(t, dtCalc, cMob);
    
    % Mass balance and transfer correctness check
    mSoluteOutStep = mSoluteOutStep + mSoluteOut;
    mSoluteOutCum = mSoluteOutCum + mSoluteOut;
    
    % Inclrease time
    t = t + dtCalc;
    
%     % Mass balance and transfer correctness check: expected theta  .* ModelDim.mobileFraction
%     thetaExp = thetaN - ...
%         (t - tPrev) * (qIn(2:end) - qIn(1:end-1)) ./ ModelDim.dzin;
%     checkExp = thetaNewN - thetaExp;
%     if (any(~RealEq(checkExp, 0, SimulationPar.EPSILON)))
%         warning('ResultValidation:MassBalance', 't = %5.3f: Fluid mass balance error: %5.3e', ...
%             t - dtCalc, max(abs(checkExp)));
%     end

    mRemainingMobile = MarkerData.dv' * MarkerData.c;
    mRemainingImmobile = (-ModelDim.dzin .* ImmobilePhase.theta)' * ImmobilePhase.c;
    mRemaining = mRemainingMobile + mRemainingImmobile;
    balance = abs(mInitial - mSoluteOutCum - mRemaining);
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
        [MarkerData, cMobOutArr(:, :, iTime)] = MarkerData.NodalConcentrations();
        cImmobOutArr(:, :, iTime) = ImmobilePhase.c;
        qOut(iTime) = vOutStep;
        mOut(:, iTime) = mSoluteOutStep;
        timeOutVec(iTime) = t;
        
        % Reset flag
        isMidStep = false;
    end
end
toc

%% Checking mass balance
fprintf('\nMass balance check:\n');
massRemainingMobile = MarkerData.dv' * MarkerData.c;
massRemainingImmobile = (-ModelDim.dzin .* ImmobilePhase.theta)' * ImmobilePhase.c;
massRemaining = massRemainingMobile + massRemainingImmobile;
balance = abs(mInitial - mSoluteOutCum - massRemaining);
for soluteIdx = 1:nSolutes
    fprintf('Solute %d: ', soluteIdx);
    fprintf('\tInitial mass: %f, flushed %f, remaining %f\n', ...
        mInitial(soluteIdx), mSoluteOutCum(soluteIdx), massRemaining(soluteIdx));
    fprintf('\t\tBalance is: %e\n', balance(soluteIdx));
    if (abs(balance(soluteIdx)) > SimulationPar.EPSILON)
        warning('ResultValidation:MassBalance', 'Final check: solute mass balance error: %5.3e', ...
            abs(balance(soluteIdx)));
    end
end
    
%% Plots
selectedNodes = [5:5:20];
DisplayPlotTwoPhases(timeOutVec, cMobOutArr, cImmobOutArr, ModelDim, SoilPar, inFlow, selectedNodes);

% Plot outflux vs. out concentration
figH = figure(2);
set(figH, 'Position', [200, 200, 400, 300]);
[axs, qHandle, mHandle] = plotyy(timeOutVec, qOut / dt, timeOutVec, mOut(soluteIdx, :) ./ qOut);
set(qHandle, 'LineStyle', '--');
set(axs(1), 'YLim', [0 -1.1 * inFlow]);
set(axs(1), 'YTick', [0:0.5:1] * -inFlow);
set(axs(2), 'YLim', [0 1]);
legend('Discharge', 'Concentration', 'Location', 'NorthEast');