classdef MarkerDataTtCl
    % TT - stands for travel time
    
	properties (Access = public)
        % Maximum volume of one marker particle
        dvMax = 1e-3;
        % Numerical diffusion coefficient
        dCoeff = 1;
        % Total number of marker particles in system (is variable)
        nTotal;
        % Number of solutes simulated
    	nSolutes;
        % Positions of markers in 1-dimensional system (dim: nTotal x 1)
        z;
        % Volumes of markers (dim: nTotal x 1)
        dv;
        % Concentrations of solutes in markers (dim: nTotal x nSolutes)
        c;
        % Velocity of markers (only in stochastic mode)
        q;
        % Diffusion coefficients of solutes (dim: 1 x nSolutes)
        d;
        % Indices of nodes to which markers belong
        node;
        % Concentrations of solutes in cells (dim: nNodes x nSolutes)
        cNodes;
        % Moisture contents in cells (dim: nNodes x 1)
        thetaNodes;
        % Fraction of mobile volume in total volume of a cell
        mobileFraction;
    end
    
    properties (Access = private)
        % Flag showing if the model is stochastic or deterministic
        isStochastic;
        % Parameters of lognormal distribution for current column
        lognormal;
        % Flag showing if markers are sorted by their location
        isSorted;
        % Flag showing if nodal concentrations have been computed recently
        hasNodalConcentrationsComputed;
        % Flag showing if nodal moisture contents have been computed recently
        hasNodalThetasComputed;
        % Numerical tolerance
        EPSILON;
        % Structures
        ModelDim;
        SoilPar;
        SimulationPar;
    end
    
    methods (Access = public)
        %% Constructor
        function self = MarkerDataTtCl(thetaN, nSolutes, ModelDim, SoilPar, SimulationPar, InitialC)
            if isfield(SimulationPar, 'isStochastic')
                self.isStochastic = SimulationPar.isStochastic;
            else
                self.isStochastic = false;
            end
            
            % Define diffusion coefficient
            self.d = SoilPar.d;
            
            % Structures
            self.ModelDim = ModelDim;
            self.SoilPar = SoilPar;
            self.SimulationPar = SimulationPar;
            
            % Tolerance
            self.EPSILON = SimulationPar.EPSILON;
            
            % Fraction of mobile volume
            if isfield(ModelDim, 'mobileFraction')
                self.mobileFraction = ModelDim.mobileFraction;
            else
                self.mobileFraction = ones(ModelDim.znn, 1);
            end
            
            % Calculate amounts of fluid in each node
            vN = thetaN .* abs(ModelDim.dzin) .* self.mobileFraction;
            
            % Calculate number of markers needed for each node
            nMarkPerN = ceil(vN / self.dvMax); % 10 * ones(ModelDim.znn, 1); % 
            
            % Initialize some numbers
            self.nTotal = sum(nMarkPerN);
            self.nSolutes = nSolutes;

            % Initialize information about markers
            self.z = zeros(self.nTotal, 1);
            self.dv = zeros(self.nTotal, 1);
            self.node = zeros(self.nTotal, 1);
            
            iPos = 1;
            % We process node by node because properties are different
            for iNode = 1:ModelDim.znn
                % Distribute markers over the node
                self.z(iPos:iPos + nMarkPerN(iNode) - 1) = ModelDim.zin(iNode) + ...
                    ModelDim.dzin(iNode) * ((1:nMarkPerN(iNode))' - 0.5) / nMarkPerN(iNode);
                % Distribute volumes of markers proportionally to moisture content
                self.dv(iPos:iPos + nMarkPerN(iNode) - 1) = ...
                    self.DistributeVolumes(vN(iNode), self.z(iPos:iPos + nMarkPerN(iNode) - 1));
                % Set indices of nodes to which markers belong
                self.node(iPos:iPos + nMarkPerN(iNode) - 1) = iNode;
                % Update index of first marker for the next node
                iPos = iPos + nMarkPerN(iNode);
            end
            
            if self.isStochastic
                self.lognormal.mu = 3;
                self.lognormal.sigma = 0.6;
                self.q = zeros(self.nTotal, 1); % self.InitializeMarkersVelocity(self.nTotal);
            end
            
            % Define initial concentration for markers
            self.c = InitialC(self.z, nSolutes);

            % Check inputs
            if (size(SoilPar.d, 2) ~= nSolutes)
                error('ParamCheck:NumSolutes', ...
                    'Specified diffusion coefficients nSolutes don''t correspond to nSolutes.');
            end
            
            % Set flags
            self.isSorted;
            self.hasNodalConcentrationsComputed = false;
            self.hasNodalThetasComputed = false;
        end
        
        %% Advect markers
        function [self, deltaT, vOut, mSoluteOut] = Advect(...
                self, t, deltaT, qIn, thetaN, BoundaryPar)
            nZn = self.ModelDim.znn;
            nZin = self.ModelDim.znin;

            thetaN = thetaN .* self.mobileFraction;
            thetaIn = InterNodalValues(self, thetaN);
            
            if ~RealGt(qIn(nZin), 0, self.EPSILON)
                % Create one extra marker with zero volume and mass at the bottom. This
                % prevents crashes when flux at the bottom is too low, and no markers leave the
                % system, but volume is expected to flow outside.
                self.nTotal = self.nTotal + 1;
                self.z      = cat(1, self.z, self.ModelDim.zin(nZin));
                self.q      = cat(1, self.q, 0);
                self.dv     = cat(1, self.dv, 0);
                self.c      = cat(1, self.c, zeros(1, self.nSolutes));
                self.node   = cat(1, self.node, nZn);
            end

            
            if ~self.isStochastic
                % Velocity of particles
                qVelIn = qIn ./ thetaIn;
                
                % Calculate maximum time step
                dzIn = -cat(1, self.ModelDim.dzin(1), self.ModelDim.dzin);
                deltaT = min(deltaT, 0.95 * min(abs(dzIn ./ qVelIn)));
                
                % Diffuse solutes
                self = Diffuse(self, t, deltaT);
                
                % Calculate flux over time interval
                qIn = qIn * deltaT;
                qVelIn = qVelIn * deltaT;
                
                % Coordinates of those points from which particles start switching to next nodes:
                zSwitch = self.ModelDim.zin - qVelIn;
                
                % Velocities of particles
                qMark = self.MarkerValues(zSwitch, qVelIn, 'linear');
                
                % Advect (dzMark = qMark)
                zNext = self.z + qMark;
                
                % Updated locations of markers
                nodeNext = self.node;
                
                % Particles that pass to next celdeltaTls:
                for iNode = 0:nZn
                    % Check direction of flow (downwards - negative, upwards - positive)
                    fluxDirection = sign(qIn(iNode + 1));
                    
                    if (fluxDirection <= 0)
                        % Flux downwards
                        % Beginning of interval where particles pass to the next node
                        isAfterSwitchPoint = RealLt(zNext, self.ModelDim.zin(iNode + 1), 0);
                        % Markers in the node from which water flows
                        isSourceNode = (self.node == iNode);
                    elseif (iNode == nZn)
                        % Leave the loop if water flows upwards at bottom
                        break
                    else
                        % Flux upwards
                        % End of interval where particles pass to current node from the next one
                        isAfterSwitchPoint = RealGe(zNext, self.ModelDim.zin(iNode + 1), 0);
                        % Markers in the node from which water flows
                        isSourceNode = (self.node == iNode + 1);
                    end
                    
                    % Indices of markers staying in current cell and moving further
                    doLeaveThisNode = isAfterSwitchPoint & isSourceNode;
                    doStayInThisNode = isSourceNode & (~doLeaveThisNode);
                    
                    % Updated nodes of markers that move forward / backward
                    nodeNext(doLeaveThisNode) = self.node(doLeaveThisNode) - fluxDirection;
                    
                    % Redistrubute fluid to keep mass of fluid flowing through internode correct
                    % Compute difference between desired and actual flux
                    diffV = fluxDirection * (qIn(iNode + 1) - ...
                        fluxDirection * sum(self.dv(doLeaveThisNode)));
                    
                    % Exchange fluid
                    if ~RealEq(diffV, 0, self.EPSILON)
                        if any(doLeaveThisNode)
                            if fluxDirection <= 0
                                iMarkSwitch = find(doLeaveThisNode, true, 'first');
                            else
                                iMarkSwitch = find(doLeaveThisNode, true, 'last');
                            end
                        else
                            if fluxDirection <= 0
                                iMarkSwitch = find(doStayInThisNode, true, 'last') + 1;
                            else
                                iMarkSwitch = find(doStayInThisNode, true, 'first') - 1;
                            end
                        end
                        
                        iMarkStay = iMarkSwitch + fluxDirection;
                        
                        self = self.LocalMassExchange(t, iMarkStay, iMarkSwitch, fluxDirection, diffV);
                    end
                end
            else
                % Determine time step
                deltaT = min(deltaT, min(abs((self.ModelDim.dzin(self.node)) ./ self.q)));
                
                % Diffuse solutes
                self = Diffuse(self, t, deltaT);
                
                % Calculate flux over time interval
                qIn = qIn * deltaT;
                
                zIn = [Inf; self.ModelDim.zin; -Inf];
                
                % Advect (dzMark = qMark)
                zNext = self.z + deltaT * self.q;
                
                % Update current node index
                dzNextNode = zIn(self.node + 2) - self.z;
                dzPrevNode = zIn(self.node + 1) - self.z;
                iNodeInc = dzNextNode >= self.q;
                iNodeDec = dzPrevNode < self.q;
                self.node(iNodeInc) = self.node(iNodeInc) + 1;
                self.node(iNodeDec) = self.node(iNodeDec) - 1;
            end

            % Update positions of markers
            self.z = zNext;
            
            % Inject new particles at the top
            if RealLt(qIn(1), 0, self.EPSILON)
                self = self.InjectFluid(t, 1, qIn, BoundaryPar.cTop, thetaIn);
            end
            % ... and bottom
            if RealGt(qIn(nZin), 0, self.EPSILON)
                self = self.InjectFluid(t, nZin, qIn, zeros(1, self.nSolutes), thetaIn);
            end
            
            % Consistency check
            self.CheckMarkersVolume(t);
            
            % Volume of fluid that leave system
            vOut = 0;
            % Masses of solutes leaving system
            mSoluteOut = zeros(1, self.nSolutes);
            
            % Remove markers that left system and calculate amounts of fluid and solutes that
            % leave the system
            isOut = (self.node > nZn);
            if any(isOut)
                % Calculate discharge rate, masses of solutes
                vOut = sum(self.dv(isOut));
                mSoluteOut(:) = self.dv(isOut)' * self.c(isOut, :);
                % Remove leaving markers
                self.z(isOut) = [];
                self.dv(isOut) = [];
                self.c(isOut, :) = [];
                self.q(isOut, :) = [];
                self.node(isOut, :) = [];
                self.nTotal = numel(self.z);
            end
            
            % Reset flags
            self.hasNodalConcentrationsComputed = false;
            self.hasNodalThetasComputed = false;
            
            % Check correctness
            if ~self.isStochastic
                self.CheckMoistureContentPerCell(t);
            end
            self.CheckConcentrations(t);
        end
        
        %% Inject fluid
        function self = InjectFluid(self, t, nodeInj, qIn, cBound, thetaIn)
            % Check direction of flow (downwards - negative, upwards - positive)
            fluxDirection = sign(qIn(nodeInj));
            % Volume of water injected
            vInj = qIn(nodeInj);
            % Point of entry
            zInj = self.ModelDim.zin(nodeInj);
            % Number of particles injected
            nMarkInj = ceil(abs(vInj) / self.dvMax);
            % Update total number of markers
            self.nTotal = self.nTotal + nMarkInj;
            
            % Distribute injected particles uniformly over the volume outside the system
            if fluxDirection < 0
                distr = (1:nMarkInj)' - 0.5;
            else
                distr = (nMarkInj:-1:1)' - 0.5;
            end
            zMarkInj = zInj - fluxDirection * distr * vInj / (thetaIn(nodeInj) * nMarkInj);
            % Calculate volumes of injected markers
            dvMarkInj = self.DistributeVolumes(abs(vInj), zMarkInj);
            % Generate random velocities
            if self.isStochastic
                qMarkInj = self.InitializeMarkersVelocity(nMarkInj);
            end
            % Append injected markers to existing arrays
            if fluxDirection > 0
                % If upwards flow - append markers at the end
                self.z = cat(1, self.z, zMarkInj);
                self.dv = cat(1, self.dv, dvMarkInj);
                self.c = cat(1, self.c, repmat(cBound, [nMarkInj, 1]));
                self.node = cat(1, self.node, nodeInj * ones(nMarkInj, 1));
                if self.isStochastic
                    self.q = cat(1, self.q, qMarkInj);
                end
            else
                % If downwards flow - append markers at the beginning
                self.z = cat(1, zMarkInj, self.z);
                self.dv = cat(1, dvMarkInj, self.dv);
                self.c = cat(1, repmat(cBound, [nMarkInj, 1]), self.c);
                self.node = cat(1, (nodeInj + 1) * ones(nMarkInj, 1), self.node);
                if self.isStochastic
                    self.q = cat(1, qMarkInj, self.q);
                end
            end
            
            % Reset flags
            self.isSorted = false;
            self.hasNodalConcentrationsComputed = false;
            self.hasNodalThetasComputed = false;
            
            if ~self.isStochastic
                % Reorder markers by coordinate
                self = self.SortMark();
                
                % Check correctness
                self.CheckMoistureContentPerCell(t);
            end
            
            self.CheckConcentrations(t);
        end
        
        %% Compute diffusion (as in Gerya's Marker-in-Cell method)
        function self = Diffuse(self, t, deltaT)
            % Mass balance checks:
            %   sum(thetaNodesX .* cNodesX .* abs(self.ModelDim.dzin))
            %   sum(self.dv .* self.c)
            
            [self, cNodesX] = NodalConcentrations(self);
            [self, thetaNodesX] = NodalThetas(self);

            tRangeOde = [0, deltaT];
            optionsOde = odeset('RelTol', 1e-5, 'AbsTol', 1e-3);
            
            cNext = zeros(self.ModelDim.znn, self.nSolutes);
            
            for soluteIdx = 1:self.nSolutes
                cNodesSol = cNodesX(:, soluteIdx);
                [~, cNodesDiff] = ode15s(...
                    @(t, cX) NettoFluxConcNodes(cX, ...
                        self.d(soluteIdx), ...
                        thetaNodesX, ...
                        self.ModelDim), ...
                        tRangeOde, ...
                        cNodesSol, ...
                        optionsOde);
                cNext(:, soluteIdx) = cNodesDiff(end, :)';
            end
           
             self.c = self.MarkerValues(self.ModelDim.zn, cNext, 'current node');
%             self = ApplySubgridDiffusion(self, t, deltaT, cNext);
            
            % Set flag
            self.hasNodalConcentrationsComputed = true;
            
            % Correctness check (set STRICT_CHECK to false - request to adjust concentrations)
            % that are outside [0; 1] interval
            self = self.CheckConcentrations(t, false);
        end
        
        %% Compute subgrid diffusion term (as in Gerya's Marker-in-Cell method)
        function self = ApplySubgridDiffusion(self, t, deltaT, cNext)
            [self, thetaNodesX] = NodalThetas(self); 
            
            % Internodal intervals for markers
            dzMark = self.MarkerValues(self.ModelDim.zn, self.ModelDim.dzin, 'current node');

            % Compute characteristic local diffusion time scale
            dtDiff = 1 ./ (2 ./ dzMark.^2 * self.d);
            
            % Concentrations at markers interpolated from updated concentrations
            cInterpMark = self.MarkerValues(self.ModelDim.zn, cNext, 'current node');
            
            % Changes in subgrid concentrations of markers
            dcSubgridMark =  ...
                (cInterpMark - self.c) .* ...
                (1 - exp(-self.dCoeff .* deltaT ./ dtDiff));
            % Changes in subgrid masses of markers
            dmSubgridMark = dcSubgridMark .* repmat(self.dv, [1, self.nSolutes]);
            % Sum up total mass change due to subgrid diffusion
            dmSubgridNode = ...
                ComputeNodalValues(self.z, dmSubgridMark, self.ModelDim, @sum, self.node);
            % Compute total needed mass change for nodes
            mNode = abs(self.ModelDim.dzin) .* thetaNodesX;
            % ... and concentrations
            dcSubgridNode = dmSubgridNode ./ repmat(mNode, [1, self.nSolutes]);
            % Based on the above compute remaining parts of concentrations for nodes
            dcRemainingNode = ...
                (cNext - self.cNodes) - dcSubgridNode;
            % ... and markers
            dcRemainingMark = ...
                self.MarkerValues(self.ModelDim.zn, dcRemainingNode, 'current node');
            % Total change in concentrations of markers
            dcMark = dcSubgridMark + dcRemainingMark;
            
            % Save calculated concentrations to object properties
            self.c = self.c + dcMark;
            self.cNodes = cNext;
        end
        
        %% Compute nodal values
        % Concentrations
        function [self, cNodes] = NodalConcentrations(self)
            if ~self.hasNodalConcentrationsComputed
                % Compute masses of solutes per each particle
                mMark = repmat(self.dv, [1, self.nSolutes]) .* self.c;
                % Compute masses of solutes per each node
                mNode = ComputeNodalValues(self.z, mMark, self.ModelDim, @sum, self.node);
                % Total volumes of solutes in nodes
                dvNode = ComputeNodalValues(self.z, self.dv, self.ModelDim, @sum, self.node);
                % Result
                cNodes = mNode ./ repmat(dvNode, [1, self.nSolutes]);
                
                % Save results to this object
                self.cNodes = cNodes;
                self.hasNodalConcentrationsComputed = true;
            end
            
            cNodes = self.cNodes;
        end
        
        % Moisture contents
        function [self, thetaRes] = NodalThetas(self)
            if ~self.hasNodalThetasComputed 
                [vN, ~] = ComputeNodalValues(self.z, self.dv, self.ModelDim, @sum, self.node);
                self.thetaNodes = - vN ./ self.ModelDim.dzin;
                self.hasNodalThetasComputed = true;
            end
            
            thetaRes = self.thetaNodes;
        end
        
        %% Set nodal values
        % Concentrations
        function self = SetNodalConcentrations(self, cNodes)
            self.cNodes = cNodes;
            self.c = self.MarkerValues(self.ModelDim.zn, cNodes, 'current node');
            
            % Update flag
            self.hasNodalConcentrationsComputed = true;
        end
    end
    
    %% Private methods
    methods (Access = private)
        % Sort particles by z coordinate
        function self = SortMark(self)
            if (~self.isSorted)
                [self.z, sIdx] = sort(self.z, 'descend');
                self.dv = self.dv(sIdx);
                self.c = self.c(sIdx, :);
                self.node = self.node(sIdx);
                self.isSorted = true;
            end
        end
        
        % Compute values at markers' locations based on nodal data
        function vMark = MarkerValues(self, zNode, vNode, INTERPOLATION_METHOD)
            if nargin < 4
                INTERPOLATION_METHOD = 'linear';
            end
            if strcmp(INTERPOLATION_METHOD, 'current node')
                vMark = vNode(self.node, :);
            else
                vMark = interp1(zNode, vNode, self.z, INTERPOLATION_METHOD, 'extrap');
            end
        end
        
        % Perform local mass exchange between markers
        function self = LocalMassExchange(self, t, iMarkStay, iMarkSwitch, fluxDirection, diff)
            %  Compute corresponding masses of solutes exchanged and exchange mass
            if (diff > 0)
                % Less fluid flows to next cell than calculated. We remove some fluid from
                % the particles that stayed and distribute it over the particle in the next cell
                iMarkDonor = iMarkStay;
                iMarkAcceptor = iMarkSwitch;
            else
                % More fluid flows than expected. We remove some fluid from the markers that
                % switched to the next node and attach it to previous particles
                iMarkDonor = iMarkSwitch;
                iMarkAcceptor = iMarkStay;
            end

            % This flag defines from which direction donors must be added ("plus" - downwards,
            % "minus" - upwards
            appendDonor = sign(diff) * fluxDirection;
            
            % Perform matterial exchange
            [self, dmTotal] = self.DetachVolume(iMarkDonor, abs(diff), appendDonor);
            self = self.AttachVolume(iMarkAcceptor, abs(diff), dmTotal);
        end
        
        % Detach total volume specified by dvTotal from markers defined by logical indexing with 
        % isMarkDonor. Returns masses of solutes removed from the system. System is updated.
        function [self, dmTotal] = DetachVolume(self, iMarkDonor, dvTotal, appendDonor)
            nMarkDonors = 1;
            vMarkDonors = self.dv(iMarkDonor);
            iLastDonor = iMarkDonor;
            % Add more donor markers if current donor marker doesn't contain enough volume
            while (vMarkDonors < dvTotal)
                iLastDonor = iLastDonor + appendDonor;
                vMarkDonors = vMarkDonors + self.dv(iLastDonor);
                nMarkDonors = nMarkDonors + 1;
            end
            iMarkDonor = iMarkDonor:appendDonor:iLastDonor;
            dvMark = dvTotal * self.dv(iMarkDonor) / vMarkDonors;
            self.dv(iMarkDonor) = self.dv(iMarkDonor) - dvMark;
            dmTotal = sum(dvMark .* self.c(iMarkDonor));
        end
        
        % Distribute volume (dvTotal) containing total masses (mTotal) of solutes over markers
        % defined with logical indices isMarkAcceptor.
        function self = AttachVolume(self, iMarkAcceptor, dvTotal, dmTotal)
            nMarkAcceptors = 1;
            dvMark = dvTotal / nMarkAcceptors;
            mMark = self.dv(iMarkAcceptor) .* self.c(iMarkAcceptor) + dmTotal / nMarkAcceptors;
            self.dv(iMarkAcceptor) = self.dv(iMarkAcceptor) + dvMark;
            self.c(iMarkAcceptor) = mMark ./ self.dv(iMarkAcceptor);
        end
        
        % Compute internodal values
        function valIn = InterNodalValues(self, valN)
            nZn = self.ModelDim.znn;
            nZin = self.ModelDim.znin;
            
            valIn = zeros(nZin, 1);
            valIn(2:nZin-1) = (valN(1:nZn-1) + valN(2:nZn)) / 2;
            valIn(1) = valN(1) + (valN(1) - valIn(2));
            valIn(nZin) = valN(nZn) + (valN(nZn) - valIn(nZin-1));
        end            
        
        % Initially distribute volumes of markers 
        %  correspondingly to moisture content for one node with given boundaries and moisture
        %  content at boundaries
        function dvIni = DistributeVolumes(self, vN, zMark)
            nMark = numel(zMark);
            dvFraction = ones(nMark, 1) / nMark;
            dvIni = vN * dvFraction;
        end
        
        % Initialize velocities of markers
        function qMark = InitializeMarkersVelocity(self, nMark)
            cdfPos = ((1:nMark)' - 0.5) ./ nMark;
            qMark = -abs(self.ModelDim.zin(self.ModelDim.znin) - self.ModelDim.zin(1)) ./ ...
                logninv(cdfPos, self.lognormal.mu, self.lognormal.sigma);
        end
        
        %% Some correctness checking procedures
        function CheckMoistureContentPerCell(self, t)
            [self, thetaNodesX] = NodalThetas(self);
            
            % Check if any cell contains more fluid than its capacity
            if any(RealGt(thetaNodesX, self.SoilPar.thetaS * self.mobileFraction, self.EPSILON))
                error('RuntimeCheck:ExceedTheta', ...
                    't = %5.3f: Moisture content is too high.', t);
            end
            
            % Check if any cell contains more fluid than its capacity
            if any(RealLt(thetaNodesX, self.SoilPar.thetaR * self.mobileFraction, self.EPSILON))
                error('RuntimeCheck:ExceedTheta', ...
                    't = %5.3f: Moisture content is too low.', t);
            end
        end
        
        function CheckMarkersVolume(self, t)
            if any(RealLt(self.dv, 0, self.EPSILON))
                error('RuntimeCheck:ExceedMarkersVolume', ...
                    't = %5.3f: Negative volume of marker(s).', t);
            end
        end
        
        function self = CheckConcentrations(self, t, STRICT_CHECK)
            if nargin < 3
                STRICT_CHECK = true;
            end
            
            % Avoid errors when negative concentrations on markers with zero volume
            isZeroVolume = RealEq(self.dv, 0, self.EPSILON);
            self.c(isZeroVolume) = 0;
            
            % Check if concentrations exceed the range [0; 1]
            isCBelowZero = RealLt(self.c(:), 0, self.EPSILON);
            isCAboveOne = RealGt(self.c(:), 1, self.EPSILON);
            
            if STRICT_CHECK
                % If strict check is requested, throw error, when concentration exceeds [0; 1]
                if any(isCBelowZero)
                    error('RuntimeCheck:ExceedConcentration', ...
                        't = %5.3f: Concentration is too low', t);
                end
                if any(isCAboveOne)
                    error('RuntimeCheck:ExceedConcentration', ...
                        't = %5.3f: Concentration is too high', t);
                end
            end
        end
        
    end                     % Private methods
end