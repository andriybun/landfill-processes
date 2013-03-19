%% TODO:    For zero and positive velocities problem exists
%%

classdef MarkerDataCl
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
        function self = MarkerDataCl(thetaN, nSolutes, ModelDim, SoilPar, SimulationPar, InitialC)
            % Fraction of mobile volume
            if isfield(ModelDim, 'mobileFraction')
                self.mobileFraction = ModelDim.mobileFraction;
            else
                self.mobileFraction = ones(ModelDim.znn, 1);
            end
            
            % Calculate amounts of fluid in each node
            vN = thetaN .* abs(ModelDim.dzin) .* self.mobileFraction;
            % Calculate number of markers needed for each node
            nMarkPerN = ceil(vN / self.dvMax);
            
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
                % Volume of each marker in node is equal to volume of fluid in node divided by
                % number of markers in that node
                self.dv(iPos:iPos + nMarkPerN(iNode) - 1) = vN(iNode) / nMarkPerN(iNode);
                % Distribute markers uniformly over node
                self.z(iPos:iPos + nMarkPerN(iNode) - 1) = ModelDim.zin(iNode) + ...
                    ModelDim.dzin(iNode) * ((1:nMarkPerN(iNode))' - 0.5) / nMarkPerN(iNode);
                % Set indices of nodes to which markers belong
                self.node(iPos:iPos + nMarkPerN(iNode) - 1) = iNode;
                % Update index of first marker for the next node
                iPos = iPos + nMarkPerN(iNode);
            end
            
            % Define initial concentration for markers
            self.c = InitialC(self.z, nSolutes);

            % Check inputs
            if (size(SoilPar.d, 2) ~= nSolutes)
                error('ParamCheck:NumSolutes', ...
                    'Specified diffusion coefficients nSolutes don''t correspond to nSolutes.');
            end
            
            % Define diffusion coefficient
            self.d = SoilPar.d;
            
            % Structures
            self.ModelDim = ModelDim;
            self.SoilPar = SoilPar;
            self.SimulationPar = SimulationPar;
            
            % Tolerance
            self.EPSILON = SimulationPar.EPSILON;
            
            % Set flags
            self.isSorted;
            self.hasNodalConcentrationsComputed = false;
            self.hasNodalThetasComputed = false;
        end
        
        %% Advect markers
        function [self, deltaT, vOut, mSoluteOut] = Advect(...
                self, t, deltaT, qIn, thetaN, BoundaryPar)
            thetaN = thetaN .* self.mobileFraction;
            thetaIn = zeros(self.ModelDim.znin, 1);
            thetaIn(2:self.ModelDim.znin-1) = ...
                (thetaN(1:self.ModelDim.znn-1) + thetaN(2:self.ModelDim.znn)) / 2;
            thetaIn(1) = thetaN(1);
            thetaIn(self.ModelDim.znin) = thetaN(self.ModelDim.znn);
                        
            % Velocity of particles
            qVelIn = qIn ./ thetaIn; %cat(1, thetaN(1), thetaN);
            
            % Calculate maximum time step
            dzIn = -cat(1, self.ModelDim.dzin(1), self.ModelDim.dzin);
            deltaT = min(deltaT, 0.95 * min(abs(dzIn ./ qVelIn)));

            % Calculate flux over time interval
            qIn = qIn * deltaT;
            qVelIn = qVelIn * deltaT;
            
            % We advect particles within ranges [zIn - q / Theta; zIn] for each internode with
            % negative velocity q and within ranges [zIn; zIn - q / Theta] for internodes with
            % positive velocities. Here positions (zIn - q / Theta) are calculated. Particles 
            % between those points are advected with linearly interpolated velocities.
            % Coordinates of those points where particles start switching to next nodes:
            zSwitch = zeros(self.ModelDim.znin, 1);
            zSwitch(1) = self.ModelDim.zin(1);
            zSwitch(2:self.ModelDim.znin) = ...
                self.ModelDim.zin(2:self.ModelDim.znin) - qVelIn(2:self.ModelDim.znin);
            zInterv = cat(1, zSwitch(2:end), self.ModelDim.zin);
            qInterv = cat(1, qVelIn(2:end), qVelIn);
            [zInterv, sIdx] = sort(zInterv, 'descend');
            qInterv = qInterv(sIdx);
            
            % Add particles with zero volume at the top and bottom of column. This is to avoid
            % errors while no particles stay in the first node or no particles leave the system
            % during time step
            self.z = cat(1, self.ModelDim.zin(1), self.z, self.ModelDim.zin(self.ModelDim.znin));
            self.dv = cat(1, 0, self.dv, 0);
            self.c = cat(1, zeros(1, self.nSolutes), self.c, zeros(1, self.nSolutes));
            self.node =cat(1, 1, self.node, self.ModelDim.znn);
            self.nTotal = self.nTotal + 2;
            
            % Velocities of particles
%             qMark = self.MarkerValues(zSwitch, qVelIn, 'linear');
            qMark = self.MarkerValues(zInterv, qInterv, 'linear');
            
            % Advect (dzMark = qMark)
            zNext = self.z + qMark;

            % Particles that pass to next cells:
            for iNode = 1:self.ModelDim.znn
                % Check direction of flow (downwards - negative, upwards - positive)
                fluxDirection = sign(qIn(iNode + 1));
                
                if (fluxDirection <= 0)
                    % Flux downwards
                    % Beginning of interval where particles pass to the next node
                    isAfterSwitchPoint = RealLt(zNext, self.ModelDim.zin(iNode + 1), 0);
                    % End of the current node
                    isBeforeNextNode = RealGe(self.z, self.ModelDim.zin(iNode + 1), 0);
                    % Index of particle that crosses internode last (with te lowest coordinate)
                    iMarkFirstSwitch = find(isAfterSwitchPoint, true, 'first');
                elseif (iNode == self.ModelDim.znn)
                    % Leave the loop if water flows upwards at bottom
                    break
                else
                    % Flux upwards
                    % End of interval where particles pass to current node from the next one
                    isAfterSwitchPoint = RealGe(zNext, self.ModelDim.zin(iNode + 1), 0);
                    % Beginning of the next node
                    isBeforeNextNode = RealLt(self.z, self.ModelDim.zin(iNode + 1), 0);
                    % Index of particle that crosses internode last (with te lowest coordinate)
                    iMarkFirstSwitch = find(isAfterSwitchPoint, true, 'last');
                end
                
                % Indices of markers staying in current cell and moving further
                % doStayInThisNode = isAfterThisNode & (~isAfterSwitchPoint);
                doLeaveThisNode = isAfterSwitchPoint & isBeforeNextNode;
                % Index of first particle that doesn't cross internode
                iMarkLastStay = iMarkFirstSwitch + fluxDirection;
                
                % Redistrubute fluid to keep mass of fluid flowing through internode correct
                % Compute difference between desired and actual flux
                diffV = fluxDirection * (qIn(iNode + 1) - ...
                    fluxDirection * sum(self.dv(doLeaveThisNode)));
                
                % Exchange fluid
                self = self.LocalMassExchange(t, iMarkLastStay, iMarkFirstSwitch, diffV);
                
                % Update nodes of markers that move forward / backward
                self.node(doLeaveThisNode) = self.node(doLeaveThisNode) - fluxDirection;
            end
            
            self.CheckMarkersVolume(t);
            
            % Update positions of markers
            self.z = zNext;
            
            % Volume of fluid that leave system
            vOut = 0;
            % Masses of solutes leaving system
            mSoluteOut = zeros(1, self.nSolutes);
            
            % Remove markers that left system and calculate amounts of fluid and solutes that
            % leave the system
            isOut = (self.z < self.ModelDim.zin(self.ModelDim.znin));
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
            else
                % Remove extra particle at the end
                self.z(end) = [];
                self.dv(end) = [];
                self.c(end) = [];
                self.node(end) = [];
                self.nTotal = self.nTotal - 1;
            end
            
            % Remove extra particle at the beginning
            if self.dv(1) == 0
                self.z(1) = [];
                self.dv(1) = [];
                self.c(1) = [];
                self.node(1) = [];
                self.nTotal = self.nTotal - 1;
            end
            
            % Inject new particles at the top and bottom
            if RealLt(qIn(1), 0, self.EPSILON)
                self = self.InjectFluid(t, 1, qIn, BoundaryPar.cTop, thetaIn);
            end
            if RealGt(qIn(self.ModelDim.znin), 0, self.EPSILON)
                self = self.InjectFluid(t, self.ModelDim.znin, ...
                    qIn, zeros(1, self.nSolutes), thetaIn);
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
        function self = InjectFluid(self, t, nodeInj, qIn, cBound, thetaIn)
            vInj = qIn(nodeInj);
            isUpwardFlow = vInj > 0;
            zInj = self.ModelDim.zin(nodeInj);
            % Number of particles injected
            nMarkInj = ceil(abs(vInj) / self.dvMax);
            % Update total number of markers
            self.nTotal = self.nTotal + nMarkInj;
            % Update volumes of markers
            self.dv = cat(1, repmat(abs(vInj) / nMarkInj, [nMarkInj, 1]), self.dv);
            % Distribute injected particles uniformly over the volume 0:vInj
            self.z = cat(1, ...
                zInj + ((1:nMarkInj) - 0.5)' * vInj / thetaIn(nodeInj) / nMarkInj, ...
                self.z);
            % Assign boundary concentration
            self.c = cat(1, repmat(cBound, [nMarkInj, 1]), self.c);
            
            % Locate new markers in first node
            self.node = cat(1, (nodeInj - isUpwardFlow) * ones(nMarkInj, 1), self.node);
            
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
        
        %% Perform local mass exchange between two neighbouring markers
        function self = LocalMassExchange(self, t, iMarkLastStay, iMarkFirstSwitch, diff)
            % Check if the first marker also changes node. This is to prevent particles from
            % moving too far forward
            if (isempty(iMarkLastStay) || iMarkLastStay == 0)
                % If yes, we check what is the balance between markers existing in current node
                % and flux to the next node. If diff is non-zero, the time stem must be
                % reduced. Otherwise mass exchange is not needed.
                if (abs(diff) > self.EPSILON)
                    error('RuntimeCheck:WrongMassExchange', ...
                        ['t = %5.3f: Local mass exchange expected ' ...
                        'with no marker to exchagne with.'], t);
                else
                    return
                end
            end
            
            %  Compute corresponding masses of solutes exchanged and exchange mass
            if (diff > 0)
                % Less fluid flows to next cell than calculated. We remove some fluid from
                % previous particle and add it to the last particle in the next cell
                iMarkDonor = iMarkLastStay;
                iMarkAcceptor = iMarkFirstSwitch;
                dvMarkAcceptor = self.dv(iMarkFirstSwitch) + diff;
                dvMarkDonor = self.dv(iMarkLastStay) - diff;
            else
                % More fluid flows than expected. We remove some fluid from the last marker
                % that changed node and attach this to previous particle
                iMarkDonor = iMarkFirstSwitch;
                iMarkAcceptor = iMarkLastStay;
                dvMarkAcceptor = self.dv(iMarkLastStay) - diff;
                dvMarkDonor = self.dv(iMarkFirstSwitch) + diff;
            end

            % Masses of solutes to be exchanged:
            mEx = self.c(iMarkDonor, :) * abs(diff);
            
            % Mass of solutes already present in target particle
            mSolMark = self.c(iMarkAcceptor, :) * self.dv(iMarkAcceptor);
            if RealEq(dvMarkAcceptor, 0, self.EPSILON)
                self.c(iMarkAcceptor, :) = 0;
            else
                self.c(iMarkAcceptor, :) = (mSolMark + mEx) / dvMarkAcceptor;
            end
            
            % Update volumes after exchange
            self.dv(iMarkAcceptor) = dvMarkAcceptor;
            self.dv(iMarkDonor) = dvMarkDonor;
            
            %% TODO: Check for recursion
            %%
            if (dvMarkDonor < 0)
                inc = iMarkDonor - iMarkAcceptor;
                self = self.LocalMassExchange(t, iMarkDonor, iMarkDonor + inc, dvMarkDonor);
                warning('RuntimeCheck:InconsistentAdvection', ...
                    't = %5.3f: badly organized advection: negative volume of marker #%d', ...
                    t, iMarkDonor);
            end
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
            else
                % Else concentrations are adjusted in order to keep mass balance correct
                if any(isCBelowZero)
                    % Indices of nodes with wrong concentrations 
                    iNodesToAdjust = unique(self.node(isCBelowZero))';
                    
                    % Process those nodes
                    for iNode = iNodesToAdjust
                        isMarkInNode = (self.node == iNode);
                        doAdjustMarker = isMarkInNode & ~(isCBelowZero);
                        massPerMarker = self.dv(doAdjustMarker) .* self.c(doAdjustMarker);
                        totalMassToRemove = sum(self.dv(isCBelowZero) .* self.c(isCBelowZero));
                        
                        % Set zero instead of negative concentrations
                        self.c(isCBelowZero) = 0;
                        
                        % Reallocate removed mass through other markers in node
                        massAdjustment = massPerMarker ./ sum(massPerMarker) .* totalMassToRemove;
                        self.c(doAdjustMarker) = ...
                            (massPerMarker + massAdjustment) ./ self.dv(doAdjustMarker);
                    end
                end
                if any(isCAboveOne)
                    % Indices of nodes with wrong concentrations 
                    iNodesToAdjust = unique(self.node(isCAboveOne))';
                    
                    % Process those nodes
                    for iNode = iNodesToAdjust
                        isMarkInNode = (self.node == iNode);
                        doAdjustMarker = isMarkInNode & ~(isCAboveOne);
                        massPerMarker = self.dv(doAdjustMarker) .* self.c(doAdjustMarker);
                        totalMassToRemove = sum(self.dv(isCAboveOne) .* (self.c(isCAboveOne) - 1));
                        
                        % Set one instead of concentrations above 1
                        self.c(isCAboveOne) = 1;
                        
                        % Reallocate removed mass through other markers in node
                        massAdjustment = massPerMarker ./ sum(massPerMarker) .* totalMassToRemove;

                        self.c(doAdjustMarker) = ...
                            (massPerMarker + massAdjustment) ./ self.dv(doAdjustMarker);
                    end
                end
            end             % if STRICT_CHECK
        end                 % Function
        
    end                     % Private methods
end