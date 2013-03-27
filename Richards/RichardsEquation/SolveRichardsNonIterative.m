function RichardsOutput = SolveRichardsNonIterative (trange,ModelDim,SoilPar,BoundaryPar,TimerPar,InitialPar)
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
dtconv = TimerPar.dtini;

% parameter controlling the picard iteration. This value is the minimum
% value for the update in pressure head (smaller values lead to improved
% mass balances, but longer simulation times)
%convcrit = TimerPar.convcrit;

%% Initialize output matrices
TimeR = tR';
hOut = hSim(:)';
qzOut = qz(:)';
qxOut=qx(:)';
thetaOut = theta(:)';
thetaOutzIN=thetazIN(:)';
thetaOutxIN=thetaxIN(:)';
CumMBal = 0;
MassBalOut = 0;
nNotConverged = 0;
CumMBalOut = CumMBal(:)';
cum_hSim = zeros(nNz,nNx);
dttOutR = dtR';
itersucces = 0;

%We are going to solve for the time derivative of the heads
%We assume that initial conditions are steady state...
hPrime = 0.*hSim(:);
hSimnp1 = 0.*hSim;

%% Time loop
while abs(tR-trange(end)) > eps,
    %our primary variable will be updated during the iteration. Store the
    %current value for the previous time
    hSim_o = hSim;
    hPrime_o = hPrime;
    cum_hSim_o = cum_hSim;
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
    
    %Check to see if it is wise to increase dtiter
    %if dtout./dtiter > 1000, 
    %    dtiter = dtout;
    %end
    
    %max time step related to flow rate
    % 50% change in pressure head
    change = 0.9;
    
    totq=-((qz(2:nINz,allnx)-qz(1:nINz-1,allnx))./repmat(dzIN,1,nNx)) -...
        ((qx(allnz,2:nINx)-qz(allnz,1:nINx-1))./repmat(dxIN,nNz,1));
    AA=theta_o(:)./totq(:);
    dtflow2=min(abs(change.*AA));
    
    while ~converged, 
        %when solution is not accepted, retry with smaller time step
        
        dtR = min([dtflow2,dtmaxR,dtconv,dtout]);
        niter = niter + 1;
        itersucces = itersucces+1;
        
        % Approximate solution:
        hSimnp1(:) = hSim(:) + dtR.*hPrime_o;
        
        Cm=Cm_function(tR,hSimnp1,SoilPar,ModelDim);
        Km = Kmat(tR,hSimnp1,SoilPar,ModelDim,BoundaryPar);
        Yv = Yvec(tR,hSimnp1,SoilPar,ModelDim,BoundaryPar);
        
        %LHS
        Atmp = (Cm - dtR.*Km);
        btmp = dtR.*Km*hSim(:)+Yv;
        
        hPrime = Atmp\btmp;
        
        delta = dtR.*(hPrime(:)-hPrime_o(:))./2;
        
        hSimtmp = hSim(:) + dtR.*hPrime(:);
        convcrit = TimerPar.reltol.*abs(hSimtmp)+TimerPar.abstol;
        
        %isp([niter tR dtR log10(max(abs(delta)))])
        
        if (max(abs(delta)-convcrit) < 0), % || niter >= maxiter),
            converged = true;
            
            hSim(:) = hSim(:) + dtR.*(hPrime_o(:)+hPrime(:))./2;
        
            %disp([tR log10(dtR) log10(max(abs(delta)))])
            
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
            
            % Time step can be calculated from the convergence criterium
            dtconv = min(2.*convcrit./abs(hPrime-hPrime_o));
            
                        
            %%Check for output
            if abs(tR-next_tout)<eps,
                TimeR = [TimeR;tR'];
                
                disp([tR log10(dtR) log10(max(abs(delta)))])
                
                %calculate time averaged hSim over output time step
                dtoutt = TimeR(end)-TimeR(end-1);
                avg_hSim = (cum_hSim-cum_hSim_o)./dtoutt;
                avghOut = [hOut;avg_hSim(:)'];
                hOut = [hOut;hSim(:)'];
                avg_States = Waterstates2D(tR,avg_hSim,SoilPar,ModelDim,BoundaryPar);
                States = Waterstates2D(tR,hSim,SoilPar,ModelDim,BoundaryPar);
                avgqzOut = [qzOut;avg_States.qz(:)'];
                avgqxOut=[qxOut;avg_States.qx(:)'];
                avgthetaOut = [thetaOut;avg_States.theta(:)'];
                avgthetaOutzIN=[thetaOutzIN;avg_States.thetazIN(:)'];
                avgthetaOutxIN=[thetaOutxIN;avg_States.thetaxIN(:)'];
                CumMBalOut = [CumMBalOut; CumMBal(:)'];
                MassBalOut=[MassBalOut;MassBal];
                dttOutR = [dttOutR;dtR'];
                qzOut = [qzOut;States.qz(:)'];
                qxOut=[qxOut;States.qx(:)'];
                thetaOut = [thetaOut;States.theta(:)'];
                thetaOutzIN=[thetaOutzIN;States.thetazIN(:)'];
                thetaOutxIN=[thetaOutxIN;States.thetaxIN(:)'];
            end
        else
            %model has not converged, retry with a smaller timestep
            dtconv = dtR.*TimerPar.DecreaseF;
                
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