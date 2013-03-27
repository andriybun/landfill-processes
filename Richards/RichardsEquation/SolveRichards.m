function RichardsOutput = SolveRichards (trange,ModelDim,SoilPar,BoundaryPar,TimerPar,InitialPar)
%% Spatial Discretization
nNx = ModelDim.nNx;
allnx = 1:nNx;
nINx = ModelDim.nINx;
xN=ModelDim.xN;
xIN=ModelDim.xIN;
dxN = ModelDim.dxN;
dxIN = ModelDim.dxIN;

nNz = ModelDim.nNz;
allnz = 1:nNz;
nINz = ModelDim.nINz;
zN = ModelDim.zN;
zIN = ModelDim.zIN;
dzN = ModelDim.dzN;
dzIN = ModelDim.dzIN;

volN = ModelDim.volN;

tic;
%% Initialize the statevariables for the simulation
tR = trange(1); %tR contains the current time for the Richards simulation
hSim = InitialPar.hIni;
States = Waterstates2D(tR,hSim,SoilPar,ModelDim,BoundaryPar);
thetazIN = States.thetazIN;
thetaxIN = States.thetaxIN;
theta = States.theta;
Seff = States.Seff;
qx = States.qx;
qz = States.qz;
Cm = Cm_function(tR,hSim,SoilPar,ModelDim);

dtmaxR = TimerPar.dtmax;
niter = 0;
maxiter = TimerPar.maxiter;
miniter = TimerPar.miniter;
maxiterRed = TimerPar.maxiterRed;
dtR = TimerPar.dtini;
dtiter = TimerPar.dtini;

% parameter controlling the picard iteration. This value is the minimum
% value for the update in pressure head (smaller values lead to improved
% mass balances, but longer simulation times)
%convcrit = TimerPar.convcrit;

CumMBal = 0;
MassBal = 0;
cum_hSim = hSim;
cum_hSim_o = hSim;
dttOutR = dtR';
hPrime = 0.*hSim;
nNotConverged = 0;


%% Initialize output matrices
nout = length(trange);

TimeR = zeros(nout,1);
dtROut = zeros(nout,1);
hOut = zeros(nout,nNz.*nNx);
avghOut = hOut;
thetaOut = zeros(nout,nNz.*nNx);
avgthetaOut = thetaOut;
qzOut = zeros(nout,nINz.*nNx);
avgqzOut = qzOut;
qxOut=zeros(nout,nNz.*nINx);
avgqxOut = qxOut;
thetaOutzIN=zeros(nout,nINz.*nNx);
avgthetaOutzIN = thetaOutzIN;
thetaOutxIN=zeros(nout,nNz.*nINx);
avgthetaOutxIN = thetaOutxIN;
MassBalOut = zeros(nout,1);
CumMBalOut = zeros(nout,1);

%Store initial values in output matrices
outCnt = 1;
TimeR(outCnt,1)=tR;
dtROut(outCnt,1) = 0;
MassBalOut (outCnt,1) = MassBal;
CumMBalOut (outCnt,1) = CumMBal;
hOut(outCnt,:) = hSim(:)';
avghOut(outCnt,:) = hSim(:)';
thetaOut(outCnt,:) = theta(:)';
avgthetaOut(outCnt,:) = theta(:)';
thetaOutzIN(outCnt,:) = thetazIN(:)';
avgthetaOutzIN(outCnt,:) = thetazIN(:)';
thetaOutxIN(outCnt,:) = thetaxIN(:)';
avgthetaOutxIN(outCnt,:) = thetaxIN(:)';
qzOut(outCnt,:) = qz(:)';
avgqzOut(outCnt,:) = qz(:)';
qxOut(outCnt,:) = qx(:)';
avgqxOut(outCnt,:) = qx(:)';


%% Time loop
while abs(tR-trange(end)) > eps,
    %our primary variable will be updated during the iteration. Store the
    %current value for the previous time
    hSim_o = hSim;
    hPrime_o = hPrime;
    theta_o = theta;
    %Theta_o=repmat(theta_o,1,nNz);
    %initialize local iteration counter (for reducing timestep)
    niter = 0;
    converged = false;
    
    %Check timestep
    %output time
    dtoutr=(tR-trange);
    dtoutidx = find(dtoutr<0);
    next_tout = trange(dtoutidx(1));
    dtout = next_tout-tR;
    
    %max time step related to flow rate
    % 50% change in pressure head
    change = 0.2;
    
    totq=-((qz(2:nINz,allnx)-qz(1:nINz-1,allnx))./repmat(dzIN,1,nNx)) -...
        ((qx(allnz,2:nINx)-qz(allnz,1:nINx-1))./repmat(dxIN,nNz,1));
    AA=theta_o(:)./totq(:);
    dtflow2=min(abs(change.*AA));
    
    while ~converged,
        dtR = min([dtflow2,dtmaxR,dtiter,dtout]);
        niter = niter + 1;
        %itersucces = itersucces+1;
        
        Cm=Cm_function(tR,hSim,SoilPar,ModelDim);
        Km = Kmat(tR,hSim,SoilPar,ModelDim,BoundaryPar);
        Yv = Yvec(tR,hSim,SoilPar,ModelDim,BoundaryPar);
        
        %LHS
        Atmp = (Cm./dtR-Km);
        btmp = Km*hSim(:)+Yv - (theta(:)-theta_o(:))./dtR;
        
        delta = Atmp\btmp;
        
        %isp([niter tR dtR log10(max(abs(delta)))])
        
        hSim(:) = hSim(:) + delta;
        
        hPrime = (hSim-hSim_o)./dtR;
        truncerr = (hPrime-hPrime_o).*dtR./2;
        convcrit = TimerPar.reltol.*abs(hSim)+TimerPar.abstol;
        testval = abs(truncerr)-convcrit;
        
        if max(testval(:)) < 0 ||  niter >= maxiter,
            converged = true;
            tR = tR+dtR;
                       
            States = Waterstates2D(tR,hSim,SoilPar,ModelDim,BoundaryPar);
            theta = States.theta;
            thetazIN=States.thetazIN;
            thetaxIN=States.thetaxIN;
            qz = States.qz;
            qx=States.qx;
    
            %Calculate cumulative states
            cum_hSim = cum_hSim+hSim*dtR;
            
            %Calculate Cumulative MassBalance
            TotWatMass_o = -sum(theta_o(:).*volN(:));
            TotWatMass = -sum(theta(:).*volN(:));
            
            % This approach to the mass balance assumes water integrates the water flux along
            % all boundaries of our rectangular domain.
            
            NFlux = -((sum(qz(1,:).*dxIN)-sum(qz(end,:).*dxIN))./(zIN(end)-zIN(1)))...
                -((sum(qx(:,1).*dzIN)-sum(qx(:,end).*dzIN))./(xIN(end)-xIN(1)));
            % Total change of mass in model domain should be equal to total the difference
            % between the time integrated flux at the top and the bottom boundaries
            IntWATmass = (TotWatMass-TotWatMass_o);
            MassBal = IntWATmass-(NFlux.*dtR);
            CumMBal = CumMBal+MassBal;
            %disp([tR MassBal CumMBal])
            %if niter >= maxiter,
            %    dtiter = dtiter./3;
            
             % Time step can be calculated from the convergence criterium
            %dtiter = min(min(abs(truncerr)./abs(hPrime)));
            
            if (niter <= miniter),% && (itersucces > TimerPar.succesiter),
                % system converged to a solution with a minimum nr of
                % timesteps, time step can be increased to speed up
                % simulation
                dtiter = dtR.*TimerPar.IncreaseF;
                itersucces = 0;
                %disp([miniter,niter tR log10(dtR) log10(max(abs(delta)))])
                
            elseif niter >= maxiter,
                % model has not converged within maxiter iterations. Model will continue after
                % reducing time step
                nNotConverged = nNotConverged+1;
                dtiter = dtR.*TimerPar.DecreaseF;
            
            elseif niter >= maxiterRed,
                % model has converged, but required a little to many iterations
                % reduce timestep
                dtiter = dtR.*TimerPar.DecreaseF;
                
                %disp([maxiter,niter tR log10(dtR) log10(max(abs(delta)))])
            end
            
            
            %%Check for output
            if abs(tR-next_tout)<eps,
                
                disp([nNotConverged tR log10(dtR) log10(max(abs(truncerr(:))))])
            
                outCnt = outCnt+1;
                TimeR(outCnt,1) = tR;
            
                %calculate time averaged hSim over output time step
                dtoutt = TimeR(outCnt)-TimeR(outCnt-1);
                avg_hSim = (cum_hSim-cum_hSim_o)./dtoutt;
                
                cum_hSim_o = cum_hSim;
                
                hOut(outCnt,:) = hSim(:)';
                avghOut = [hOut;avg_hSim(:)'];
                
                avg_States = Waterstates2D(tR,avg_hSim,SoilPar,ModelDim,BoundaryPar);
                States = Waterstates2D(tR,hSim,SoilPar,ModelDim,BoundaryPar);
                
                qzOut(outCnt,:) = States.qz(:)';
                qxOut(outCnt,:) = States.qx(:)';
                thetaOut(outCnt,:) = States.theta(:)';
                thetaOutzIN(outCnt,:) = States.thetazIN(:)';
                thetaOutxIN(outCnt,:) = States.thetaxIN(:)';
                
                avgqzOut(outCnt,:) = avg_States.qz(:)';
                avgqxOut(outCnt,:) = avg_States.qx(:)';
                avgthetaOut(outCnt,:) = avg_States.theta(:)';
                avgthetaOutzIN(outCnt,:) = avg_States.thetazIN(:)';
                avgthetaOutxIN(outCnt,:) = avg_States.thetaxIN(:)';
                
                
                CumMBalOut(outCnt,1) = CumMBal;
                MassBalOut = MassBal;
                dtROut(outCnt,1) = dtR;
                
            end
        end
    end
end
toc;
RichardsOutput.trange = trange;
RichardsOutput.dttOutR = dttOutR;
RichardsOutput.avghOut = avghOut;
RichardsOutput.hOut = hOut;
RichardsOutput.avgthetaOut = avgthetaOut; 
RichardsOutput.thetaOut = thetaOut;
RichardsOutput.avgthetaOutzIN = avgthetaOutzIN;
RichardsOutput.thetaOutzIN = thetaOutzIN;
RichardsOutput.avgthetaOutxIN = avgthetaOutxIN;
RichardsOutput.thetaOutxIN = thetaOutxIN;
RichardsOutput.avgqzOut = avgqzOut;
RichardsOutput.qzOut = qzOut;
RichardsOutput.avgqxOut = avgqxOut;
RichardsOutput.qxOut = qxOut;
RichardsOutput.MassBalOut = MassBalOut;
RichardsOutput.CumMBalOut = CumMBalOut;
RichardsOutput.nNotConverged = nNotConverged;

save RichardsResults RichardsOutput;