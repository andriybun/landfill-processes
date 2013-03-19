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
            
            %% TODO: initialize concentration instead of this:
            %%
            self.c = ones(self.nNodes, self.nSolutes);

            % Structures
            self.ModelDim = ModelDim;
            self.SoilPar = SoilPar;
            self.SimulationPar = SimulationPar;
            self.ModelDim.nPhases = self.nPhases;
            
            % Tolerance
            self.EPSILON = SimulationPar.EPSILON;
        end
        
        %% Compute diffusion
        function [self, cMobDiff] = Diffuse(self, t, deltaT, cMob, thetaMob)
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
                    
%                     % Replace NaN's where theta was zero
%                     isThetaZero = (xTheta == 0);
%                     cNextSol(isThetaZero) = 0;
                    
                    cDiff(doComputeNodes, :, iSolute) = cNextSol;
                    
                    % Mass balance check
                    chk = ~RealEq(sum(cSol .* xTheta, 2), sum(cNextSol .* xTheta, 2), self.EPSILON);
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
        end
    end
    
    %% Private methods
    methods (Access = private)
        
    end
    
    %% Static methods
    methods (Static)
        
    end
end