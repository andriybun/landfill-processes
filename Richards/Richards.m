function [qOut, thetaOut, hOut, time] = Richards(tRange, dtMax, ModelDim, SoilPar, BoundaryPar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Richards Equation:
% Equation form: 
%   cMatr.dh/dt = -d(Ksat.kr)/dz * (dh/dz + 1))
%
% where the dependent variable is h-pressure head
%   cMatr:   Specific mositure capacity
%   z:    spatial dimension, z direction
%   Ksat: saturated hydraulic conductivity
%   kr:   relative permeability
%   t:    temporal dimesion, time
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    allNodes = 1:ModelDim.znn;
    allInodes = 1:ModelDim.znin;

    % Model initial states;
    rhow = 1000;        % kg/m3 density of water
    g = 9.81;           % m/s^2 gravitational acceleration
    hIni(allNodes, 1) = BoundaryPar.hBot - (ModelDim.zn - BoundaryPar.zRef);

    %% Parameters for iterations
    maxIter = 45;
    minIter = 15;
    dt = dtMax;
    dtIter = dtMax;

    CONVERGENCE_EPSILON = 1e-7;

    % Initialize state variables
    t = tRange(1);
    hSim = hIni;
    states = ComputeWaterStates(t, hSim, SoilPar, ModelDim, BoundaryPar);
    theta = states.theta;
    q = states.q;
    cMatr = ComputeSpecificMoistureCapacityMatrix(hSim, SoilPar);

    % Initialize output matrices
    time = t;
    hOut = hSim';
    qOut = q';
    thetaOut = theta';
    cumMassBalance = 0;
    cumMassBalanceOut = cumMassBalance;
    dttOut = dt;
    
    % Initialize cumulative flux at current time interval. This is used to calculate time
    % averaged flux based on all sub-intervals
    qCum = zeros(ModelDim.znin, 1);
    dtCum = 0;
    
    while abs(t - tRange(end)) > eps
        % Our primary variable will be updated during the iteration. Store the
        % current value for the previous time
        hSimPrev = hSim;
        thetaPrev = theta;

        % Initialize local iteration counter (for reducing timestep)
        nIter = 0;
        hasConverged = false;

        % Check timestep
        % Output time
        dtOutR = (t - tRange);
        dtOutIdx = find(dtOutR < 0);
        tOutNext = tRange(dtOutIdx(1));
        dtout = tOutNext-t;

        % max time step related to flow rate
        % 50% change in pressure head
        change = 0.9;
        dtflow1 = min(abs((change.*cMatr*hSim.*ModelDim.dzin)./-(q(2:ModelDim.znin, 1)-q(1:ModelDim.znin-1, 1))));
        dtflow2 = min(abs((change.*thetaPrev.*ModelDim.dzin)./-(q(2:ModelDim.znin, 1)-q(1:ModelDim.znin-1, 1))));

        dt = min([dtflow2, dtMax, dtIter, dtout]);

        while ~hasConverged
            nIter = nIter + 1;
            cMatr = ComputeSpecificMoistureCapacityMatrix(hSim, SoilPar);
            kMatr = ComputeConductivityMatrix(hSim, SoilPar, ModelDim, BoundaryPar);
            yVec = ComputeY(t, hSim, SoilPar, ModelDim, BoundaryPar);

            % (kMatr + cMatr) * delta = kMatr * h + yVec - (theta - thetaPrev) ./ dt
            % Ax = b --> x = A \ b
            aMatr = (-kMatr + cMatr ./ dt);
            bVec = kMatr * hSim + yVec - (theta - thetaPrev) ./ dt;
            deltaH = aMatr \ bVec;

            hSim = hSim + deltaH;
            States = ComputeWaterStates(t, hSim, SoilPar, ModelDim, BoundaryPar);
            theta = States.theta;
            q = States.q;
            
            if max(abs(deltaH)) < CONVERGENCE_EPSILON
                hasConverged = true;
                if nIter >= maxIter
                    dtIter = dtIter ./ 3;
                elseif nIter <= minIter
                    dtIter = dtIter .* 2;
                end
                t = t + dt;
                
                % Calculate Cumulative Mass Balance
                totalWaterMassPrev = -sum(thetaPrev(:, 1) .* ModelDim.dzin);
                totalWaterMass = -sum(theta(:, 1) .* ModelDim.dzin);
                nettoFlux = -(q(1, 1) - q(end, 1));
                
                % Total change of mass in column should be equal to total the 
                % difference between the time integrated flux at the top and 
                % the bottom boundaries
                massBallance = (totalWaterMass - totalWaterMassPrev) - nettoFlux .* dt;
                cumMassBalance = cumMassBalance + abs(massBallance);

                % Accumulate fluxes and time steps
                qCum = qCum + q * dt;
                dtCum = dtCum + dt;
                
                % Check for output
                if abs(t - tOutNext) < eps
                    time = [time; t];
                    hOut = [hOut; hSim'];
                    qOut = [qOut; qCum' / dtCum];
                    
                    % Reset cumulative flux at current time interval
                    qCum = zeros(ModelDim.znin, 1);
                    dtCum = 0;
                    
                    thetaOut = [thetaOut; theta'];
                    cumMassBalanceOut = [cumMassBalanceOut; cumMassBalance];
                    dttOut = [dttOut; dt];
                end

            end
        end     % while not has converged
    end     % while - time loop
    
end     % function