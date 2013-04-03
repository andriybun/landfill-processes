% Marker Data 2D class
classdef MarkerDataCl
	properties (Access = public)
        % Maximum volume of one marker particle
        dvMax = 1e-4;
        % Numerical diffusion coefficient
        dCoeff = 1;
        % Total number of marker particles in system (is variable)
        nTotal;
        % Number of solutes simulated
    	nSolutes;
        % Positions of markers in 2-dimensional system (dim: nTotal x 1)
        z;
        x;
        % Volumes of markers (dim: nTotal x 1)
        dv;
        % Concentrations of solutes in markers (dim: nTotal x nSolutes)
        c;
        % Diffusion coefficients of solutes (dim: 1 x nSolutes)
        d;
        % Indices of nodes to which markers belong
        zNode;
        xNode;
        % Concentrations of solutes in cells (dim: nNodes x nSolutes)
        cNodes;
        % Moisture contents in cells (dim: nNodes x 1)
        thetaNodes;
        % Fraction of mobile volume in total volume of a cell
        mobileFraction;
    end
    
    properties (Access = private)
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
        
        % Some constants
        FROM_ABOVE = 0;
        FROM_BELOW = 1;
    end
    
    methods (Access = public)
        %% Constructor
        function self = MarkerDataCl(thetaN, nSolutes, ModelDim, SoilPar, SimulationPar, InitialC)
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
                self.mobileFraction = ones(ModelDim.znn, ModelDim.xnn);
            end
            
            % Calculate amounts of fluid in each node
            vN = thetaN .* (abs(ModelDim.dzin) * abs(ModelDim.dxin)) .* self.mobileFraction;
            
            % Calculate number of markers needed for each node
            nMarkPerNode = round(vN / self.dvMax);
            
            % Initialize some numbers
            self.nTotal = sum(sum(nMarkPerNode));
            self.nSolutes = nSolutes;

            % Initialize information about markers
            self.z = zeros(self.nTotal, 1);
            self.x = zeros(self.nTotal, 1);
            self.dv = zeros(self.nTotal, 1);
            self.zNode = zeros(self.nTotal, 1);
            self.xNode = zeros(self.nTotal, 1);
            
            iPos = 1;
            % We process node by node because properties are different
            for iNodeZ = 1:ModelDim.znn
                for iNodeX = 1:ModelDim.xnn
                    nMark = nMarkPerNode(iNodeZ, iNodeX);
                    % Allocate markers in the node
                    [self.z(iPos:iPos + nMark - 1), self.x(iPos:iPos + nMark - 1)] = ...
                        self.AllocateMarkers(...
                        [ModelDim.zin(iNodeZ), ModelDim.zin(iNodeZ + 1)], ...
                        [ModelDim.xin(iNodeX), ModelDim.xin(iNodeX + 1)], ...
                        nMark);
                    % Distribute volumes of markers proportionally to moisture content in cell
                    self.dv(iPos:iPos + nMark - 1) = ...
                        self.DistributeVolumes(vN(iNodeZ, iNodeX), nMark);
                    % Set indices of nodes to which markers belong
                    self.zNode(iPos:iPos + nMark - 1) = iNodeZ;
                    self.xNode(iPos:iPos + nMark - 1) = iNodeX;
                    % Update index of first marker for the next node
                    iPos = iPos + nMark;
                end
            end
            
            % Define initial concentration for markers
            self.c = InitialC(self.z, self.x, nSolutes);

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
                self, t, deltaT, qzIn, qxIn, thetaN, BoundaryPar)
            nZn = self.ModelDim.znn;
            nZin = self.ModelDim.znin;
            nXn = self.ModelDim.xnn;
            nXin = self.ModelDim.xnin;

            thetaN = thetaN .* self.mobileFraction;
            [thetaZin, thetaXin] = InterNodalValues(self, thetaN);
                        
            % Velocity of particles (vertical and horizontal components)
            qzVelIn = qzIn ./ thetaZin;
            qxVelIn = qxIn ./ thetaXin;
            
            % Calculate maximum time step
            dzIn = -cat(1, self.ModelDim.dzin(1), self.ModelDim.dzin);
            dxIn = -cat(2, self.ModelDim.dxin(1), self.ModelDim.dxin);
            deltaT = min([deltaT, ...
                0.95 * min(min(abs(repmat(dzIn, [1, nXn]) ./ qzVelIn))), ...
                0.95 * min(min(abs(repmat(dxIn, [nZn, 1]) ./ qxVelIn)))]);

            % Calculate flux over time interval
            qzIn = qzIn * deltaT;
            qxIn = qxIn * deltaT;
            qzVelIn = qzVelIn * deltaT;
            qxVelIn = qxVelIn * deltaT;
            
            % Inject new particles at the top
            self = self.InjectFluidVertical(t, 1, self.FROM_ABOVE, qzIn, ...
                BoundaryPar.cTop, thetaZin);
            % ... and bottom
            self = self.InjectFluidVertical(t, nZn, self.FROM_BELOW, qzIn, ...
                zeros(1, self.nSolutes), thetaZin);
            
            % Coordinates of those points from which particles start switching to next nodes:
            zSwitch = self.ModelDim.zin - qVelIn;
            
            % Velocities of particles
            qMark = self.MarkerValues(zSwitch, qVelIn, 'linear');
            
            % Advect (dzMark = qMark)
            zNext = self.z + qMark;

            % Updated locations of markers
            nodeNext = self.node;
            
            % Particles that pass to next cells:
            for iNode = 0:nZn
                % Check direction of flow (downwards - negative, upwards - positive)
                fluxDirection = sign(qzIn(iNode + 1));
                
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
                diffV = fluxDirection * (qzIn(iNode + 1) - ...
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

            % Update positions of markers
            self.z = zNext;
            self.node = nodeNext;
            
%             % Consistency check
%             self.CheckMarkersVolume(t);
            
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
                self.node(isOut, :) = [];
                self.nTotal = numel(self.z);
            end
            
            % Reset flags
            self.hasNodalConcentrationsComputed = false;
            self.hasNodalThetasComputed = false;
            
            % Diffuse solutes
            self = Diffuse(self, t, deltaT);
            
            % Check correctness
            self.CheckMoistureContentPerCell(t);
            self.CheckConcentrations(t);
        end
        
        %% Inject fluid
        function self = InjectFluidVertical(self, t, nodeInj, sourcePos, qzIn, cBound, thetaZin)
            % nodeInj       - index of node to inject to
            % sourcePos     - 0 - from above, 1 - from below

            % Index of inode through which injection is done
            iNodeInj = nodeInj + sourcePos;
            % Check direction of flow (positive - inside, negative - outside)
            fluxDirection = sign(qzIn(iNodeInj, :)) * sign(self.ModelDim.dzin(nodeInj)) * ...
                (2 * (sourcePos - 0.5));
            % Check for influx
            doInject = (fluxDirection < 0);
            
            % Volume of water injected
            vInj = zeros(1, self.ModelDim.xnn);
            vInj(doInject) = abs(qzIn(iNodeInj, doInject));
            
            % Point of entry
            zInj = self.ModelDim.zin(iNodeInj);
            % Number of particles injected (at least 1)
            nMarkInj = max(1, ceil(vInj / self.dvMax));
            % Update total number of markers
            self.nTotal = self.nTotal + sum(nMarkInj);
            
            % Distribute injected particles uniformly over the volume outside the system
            zMarkInj = zeros(sum(nMarkInj), 1);
            xMarkInj = zeros(sum(nMarkInj), 1);
            
            iSt = 1;
            for iX = 1:self.ModelDim.xnn
                iEnd = iSt + nMarkInj(iX) - 1;
                [zMarkInj(iSt:iEnd), xMarkInj(iSt:iEnd)] = self.AllocateMarkers(...
                    [zInj, zInj + (2 * (sourcePos - 0.5) * sign(self.ModelDim.dzin(nodeInj)) * ...
                    vInj(iX) / thetaZin(iNodeInj, iX))], ...
                    [self.ModelDim.xin(iX), self.ModelDim.xin(iX + 1)], ...
                    nMarkInj(iX));
            end
            
            % Calculate volumes of injected markers
            dvMarkInj = self.DistributeVolumes(vInj, nMarkInj);
            
            % Append injected markers to existing arrays
            % If downwards flow - append markers at the beginning
            self.z = cat(1, zMarkInj, self.z);
            self.x = cat(1, xMarkInj, self.x);
            self.dv = cat(1, dvMarkInj, self.dv);
            self.c = cat(1, repmat(cBound, [sum(nMarkInj), 1]), self.c);
            self.node = cat(1, (nodeInj + fluxDirection) * ones(nMarkInj, 1), self.node);
            
            % Reset flags
            self.isSorted = false;
            self.hasNodalConcentrationsComputed = false;
            self.hasNodalThetasComputed = false;
            
            % Reorder markers by coordinate
            self = self.SortMark();
            
            % Check correctness
            self.CheckMoistureContentPerCell(t);
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
           
            self = ApplySubgridDiffusion(self, t, deltaT, cNext);
            
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
                mNode = ComputeNodalValues(self.z, self.x, mMark, self.ModelDim, @sum, ...
                    self.zNode, self.xNode);
                % Total volumes of solutes in nodes
                dvNode = ComputeNodalValues(self.z, self.x, self.dv, self.ModelDim, @sum, ...
                    self.zNode, self.xNode);
                % Result
                cNodes = mNode ./ repmat(dvNode, [1, 1, self.nSolutes]);
                
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
    
    %% 
    methods (Access = private)
        %% Sort particles by z coordinate
        function self = SortMark(self)
            if (~self.isSorted)
                [self.z, sIdx] = sort(self.z, 'descend');
                self.dv = self.dv(sIdx);
                self.c = self.c(sIdx, :);
                self.node = self.node(sIdx);
                self.isSorted = true;
            end
        end
        
        %% Compute values at markers' locations based on nodal data
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
        
        %% Perform local mass exchange between markers
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
        
        %% Compute internodal values
        function [valZin, valXin] = InterNodalValues(self, valN)
            nZn = self.ModelDim.znn;
            nZin = self.ModelDim.znin;
            nXn = self.ModelDim.xnn;
            nXin = self.ModelDim.xnin;
            
            valZin = zeros(nZin, nXn);
            valZin(2:nZin-1, :) = (valN(1:nZn-1, :) + valN(2:nZn, :)) / 2;
            valZin(1, :) = valN(1, :) + (valN(1, :) - valZin(2, :));
            valZin(nZin, :) = valN(nZn, :) + (valN(nZn, :) - valZin(nZin-1, :));
            
            valXin = zeros(nZn, nXin);
            valXin(:, 2:nXin-1) = (valN(:, 1:nXn-1) + valN(:, 2:nXn)) / 2;
            valXin(:, 1) = valN(:, 1) + (valN(:, 1) - valXin(:, 2));
            valXin(:, nXin) = valN(:, nXn) + (valN(:, nXn) - valXin(:, nXin-1));
        end            
        
        %% Initially distribute volumes of markers 
        %  correspondingly to moisture content for one node with given boundaries and moisture
        %  content at boundaries
        function dvIni = DistributeVolumes(self, vN, nMark)
            dvFraction = ones(nMark, 1) / nMark;
            dvIni = vN * dvFraction;
        end
         
        %% Distribute markers' positions randomly according to latin hypercube sampling
        function [zPos, xPos] = AllocateMarkers(self, zRange, xRange, nMark)
            coords = lhspoint(nMark, 2);
            zPos = coords(1, :) * (zRange(2) - zRange(1)) + zRange(1);
            xPos = coords(2, :) * (xRange(2) - xRange(1)) + xRange(1);
        end
        
        %% Some correctness checking procedures
        function CheckMoistureContentPerCell(self, t)
            [self, thetaNodesX] = NodalThetas(self);
            
            % Check if any cell contains more fluid than its capacity
            if any(RealGt(thetaNodesX, self.SoilPar.thetaS, self.EPSILON))
                error('RuntimeCheck:ExceedTheta', ...
                    't = %5.3f: Moisture content is too high.', t);
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
%             else
%                 % Else concentrations are adjusted in order to keep mass balance correct
%                 if any(isCBelowZero)
%                     % Indices of nodes with wrong concentrations 
%                     iNodesToAdjust = unique(self.node(isCBelowZero))';
%                     
%                     % Process those nodes
%                     for iNode = iNodesToAdjust
%                         isMarkInNode = (self.node == iNode);
%                         doAdjustMarker = isMarkInNode & ~(isCBelowZero);
%                         massPerMarker = self.dv(doAdjustMarker) .* self.c(doAdjustMarker);
%                         totalMassToRemove = sum(self.dv(isCBelowZero) .* self.c(isCBelowZero));
%                         
%                         % Set zero instead of negative concentrations
%                         self.c(isCBelowZero) = 0;
%                         
%                         % Reallocate removed mass through other markers in node
%                         massAdjustment = totalMassToRemove / numel(massPerMarker);
%                         self.c(doAdjustMarker) = ...
%                             (massPerMarker + massAdjustment) ./ abs(self.dv(doAdjustMarker));
%                     end
%                 end
%                 if any(isCAboveOne)
%                     % Indices of nodes with wrong concentrations 
%                     iNodesToAdjust = unique(self.node(isCAboveOne))';
%                     
%                     % Process those nodes
%                     for iNode = iNodesToAdjust
%                         isMarkInNode = (self.node == iNode);
%                         doAdjustMarker = isMarkInNode & ~(isCAboveOne);
%                         massPerMarker = self.dv(doAdjustMarker) .* self.c(doAdjustMarker);
%                         totalMassToRemove = sum(self.dv(isCAboveOne) .* (self.c(isCAboveOne) - 1));
%                         
%                         % Set one instead of concentrations above 1
%                         self.c(isCAboveOne) = 1;
%                         
%                         % Reallocate removed mass through other markers in node
%                         massAdjustment = totalMassToRemove / numel(massPerMarker);
% 
%                         self.c(doAdjustMarker) = ...
%                             (massPerMarker + massAdjustment) ./ abs(self.dv(doAdjustMarker));
%                     end
%                 end
            end             % if STRICT_CHECK
        end                 % Function
        
    end                     % Private methods
end