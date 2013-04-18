classdef ImmobilePhaseCl
    %% Properties
    properties (Access = public)
        % Total number of nodes in system
        nNodes;
        % Number of solutes simulated
    	nSolutes;
        % Number of phases simulated
        nPhases = 2;
        % Nodal theta in immobile phase
        theta;
        % Nodal concentrations in immobile phase
        c;
        % Fraction of mobile volume in total volume of a cell
        mobileFraction;
    end
    
    properties (Access = private)
        % Densities of solutes
        rho;
        % Initial masses of solutes in solid state
        mSolidIni;
        % Decay rates
        lambda;
        % Numerical tolerance
        EPSILON;
        % Structures
        ModelDim;
        SoilPar;
        SimulationPar;
    end
    
    %% Public methods
    methods (Access = public)
        %% Constructor
        function self = ImmobilePhaseCl(theta, ...
                                        nSolutes, ...
                                        ModelDim, ...
                                        SoilPar, ...
                                        SimulationPar)
            self.nNodes = ModelDim.znn;
            self.nSolutes = nSolutes;
            
            self.mobileFraction = ModelDim.mobileFraction;
            self.theta = theta .* (1 - self.mobileFraction);
            
            %% TODO: initialize concentrations and densities instead of this:
            %%
            self.c = ones(self.nNodes, self.nSolutes);
            self.rho = ones(1, self.nSolutes + 1); % nSolutes + 1 for intert materails
            soluteFractions = 1 / (self.nSolutes + 1) * ones(self.nNodes, self.nSolutes + 1);
            %% TODO: decay constant
            self.lambda = 0 * 1e-1 * ones(1, self.nSolutes); 
            %%
            
            % Structures
            self.ModelDim = ModelDim;
            self.SoilPar = SoilPar;
            self.SimulationPar = SimulationPar;
            self.ModelDim.nPhases = self.nPhases;
            
            % Tolerance
            self.EPSILON = SimulationPar.EPSILON;
            
            % Initialize masses in solid phase
            self = self.InitializeSolidPhase(soluteFractions);
        end
        
        %% Compute diffusion
        function [self, cMobDiff] = Diffuse(self, t, deltaT, cMob, thetaMob)
            % Dissolve some matter from solid to liquid
            self = DissolveSolutes(self, t, deltaT);
            
            % cMob      - concentration of solutes in mobile phase (dim: nNodes x nSolutes)
            % thetaMob  - moisture content in mobile phase(dim: nNodes x 1)

            xTheta = cat(2, thetaMob, self.theta);
            doComputeNodes = ~any(RealEq(xTheta', 0, self.EPSILON))';
            cPhases = cat(2, permute(cMob, [1, 3, 2]), permute(self.c, [1, 3, 2]));
            
            xSoilPar = self.SoilPar;
            xSoilPar.theta = xTheta(doComputeNodes, :);
            
            xModelDim = self.ModelDim;
            xModelDim.znn = sum(doComputeNodes);
            
            if xModelDim.znn == 0
                cMobDiff = cMob;
            else
                % Array for resulting concentrations
                cDiff = zeros(xModelDim.znn, self.nPhases, self.nSolutes);
                
                % Process each solute separately (in order to be able to work with ODE solver
                for iSolute = 1:self.nSolutes
                    % Here we solve equation:
                    %   d(Theta * C) / dt = kExch * (Cmob - Cimmob)
                    cSol = cPhases(doComputeNodes, :, iSolute);
                    
                    % ODE parameters and solution
                    tRangeOde = [0, deltaT];
                    options = odeset('RelTol', 1e-3, 'AbsTol', 1e-3);
                    [~, cNextSol] = ode15s(...
                        @(t, cX) NettoFluxPhaseDiffusion(t, cX, xModelDim, xSoilPar), ...
                        tRangeOde, cSol, options);
                    
                    % Now we have to translate solutiion to concentration:
                    cNextSol = reshape(cNextSol(end, :), [], self.nPhases);
                    
                    cDiff(doComputeNodes, :, iSolute) = cNextSol;
                    
                    % Mass balance check
                    chk = ~RealEq(sum(cSol .* xTheta(doComputeNodes, :), 2), ...
                        sum(cNextSol .* xTheta(doComputeNodes, :), 2), self.EPSILON);
                    iProblemNode = find(chk, 1, 'first');
                    if iProblemNode
                        warning('RuntimeCheck:MassBalanceError', ...
                            't = %5.3f: Mass balance problem in %d-th node (solute#%d\n', ...
                            t, iProblemNode, iSolute);
                    end
                    % End mass balance check
                end
                
                % Process results
                cMobDiff = squeeze(cDiff(:, 1, :));
                self.c = squeeze(cDiff(:, 2, :));
            end
            
            self.CheckConcentrations(t);
        end
        
        %% Dissolve solutes from solid state to liquid in immobile phase
        function self = DissolveSolutes(self, t, deltaT)
            dmSolid = repmat(self.lambda, [self.nNodes, 1]) .* ...
                self.mSolidIni(:, 1:self.nSolutes) .* deltaT;
            % Remove this mass from solid phase
            self.mSolidIni(:, 1:self.nSolutes) = ...
                self.mSolidIni(:, 1:self.nSolutes) - dmSolid;
            % Add it to liquid phase
            dvLiquid = dmSolid ./ repmat(self.rho(1:self.nSolutes), [self.nNodes, 1]);
            dTheta = sum(dvLiquid, 2) ./ abs(self.ModelDim.dzin);
            vNodeSol = repmat(abs(self.ModelDim.dzin) .* self.theta, [1, self.nSolutes]) .* self.c;
            vNodeSol = vNodeSol + dvLiquid;
            self.theta = self.theta + dTheta;
            self.c = vNodeSol ./ repmat(abs(self.ModelDim.dzin) .* self.theta, [1, self.nSolutes]);
        end
    end
    
    %% Private methods
    methods (Access = private)
        function self = InitializeSolidPhase(self, soluteFractions)
            % Initial masses of different compounds in solid phase
            mSolidTotal = ...
                abs(self.ModelDim.dzin) .* (1 - self.SoilPar.thetaS) ./ ...
                sum(soluteFractions ./ repmat(self.rho, [self.nNodes, 1]), 2);
            self.mSolidIni = repmat(mSolidTotal, [1, self.nSolutes + 1]) .* soluteFractions;
        end
        
        %% Some correctness checking procedures
        function CheckConcentrations(self, t)
            % Check if concentrations exceed the range [0; 1]
            isCBelowZero = RealLt(self.c(:), 0, self.EPSILON);
            isCAboveOne = RealGt(self.c(:), 1, self.EPSILON);
            
            
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
    
    %% Static methods
    methods (Static)
        
    end
end