function nF = NettoFluxPhaseDiffusion(t, c, ModelDim, SoilPar)
    nPhases = ModelDim.nPhases;
    
    % Phases:
    %   1 - mobile
    %   2 - immobile
    
    % Calculate relative flux
    cAbs = reshape(c, [], nPhases);
    nF = SoilPar.kExch .* (cAbs(:, 2) - cAbs(:, 1));
    
    % Adjust flux of concentration based on ratio between moisture (mobile and immobile) contents
    nF = nF .* min(SoilPar.theta, [], 2);
    nF = cat(2, nF ./ SoilPar.theta(:, 1), -nF ./ SoilPar.theta(:, 2));
    
    % Reshape for ODE solver
    nF = reshape(nF, [], 1);
    
end