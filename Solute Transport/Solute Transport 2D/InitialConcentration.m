function cIni = InitialConcentration(z, x, nSolutes)
    if nargin == 1
        nSolutes = 1;
    end
    cIni = ones(numel(z), nSolutes);
    % Concentration = 1 for z in [0, -1] and x <= 0.5
    cIni = repmat((z > -1) & (z <= 0) & (x <= 0.5), [1, nSolutes]);
end