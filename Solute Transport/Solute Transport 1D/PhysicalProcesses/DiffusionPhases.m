function cDiff = DiffusionPhases(t, ModelDim, SoilPar, SimulationPar)
    %
    % Function to compute diffusion between phases
    %
    
    EPSILON = SimulationPar.EPSILON;
    
    nPhases = ModelDim.nPhases;
    nSolutes = SoilPar.nSolutes;
    
    theta = SoilPar.theta;
    
    % Initialize output
    cDiff = zeros(ModelDim.znn, nPhases, nSolutes);
    
    % Process each solute separately (in order to be able to work with ODE solver
    for iSolute = 1:nSolutes
        % Equation:
        %   d(Theta * C) / dt = div(Theta * grad(C))
        % will be solved. Hence we solve it with respect to Theta * C
        cSol = theta .* SoilPar.c(:, :, iSolute);
        d = SoilPar.d(iSolute);

        % ODE parameters and solution
        tRangeOde = [t, t + SimulationPar.dt];
        options = odeset('RelTol', 1e-3, 'AbsTol', 1e-3);
        [~, cThetaNextSol] = ode15s(...
            @(t, cSol) NettoFluxPhaseDiffusion(t, cSol, ModelDim, SoilPar, d), ...
            tRangeOde, cSol, options);
        
        % Now we have to translate solutiion to concentration:
        cThetaNextSol = reshape(cThetaNextSol(end, :), [], nPhases) ./ theta;
        
        % Replace NaN's where theta was zero
        isThetaZero = (theta == 0);
        cThetaNextSol(isThetaZero) = 0;
        
        cDiff(:, :, iSolute) = cThetaNextSol;

        % Mass balance check
        chk = (abs(sum(cSol, 2) - sum(cThetaNextSol .* theta, 2)) >= EPSILON);
        iProblemNode = find(chk, 1, 'first');
        if iProblemNode
            warning('Matlab:MassBalance', ...
                'Mass balance problem in %d-th node (solute#%d0\n', iProblemNode, iSolute);
        end
        % End mass balance check
    end
end