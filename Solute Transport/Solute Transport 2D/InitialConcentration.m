function cIni = InitialConcentration(z, nSolutes)
    if nargin == 1
        nSolutes = 1;
    end
    cIni = 0 * ones(numel(z), nSolutes);
    % Concentration = 1 in interval [0, -1]
    cIni = repmat((z > -1) & (z <= 0), [1, nSolutes]);
end