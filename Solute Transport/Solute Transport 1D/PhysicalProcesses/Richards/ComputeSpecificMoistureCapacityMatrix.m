function cMatr = ComputeSpecificMoistureCapacityMatrix(h, SoilPar)
    alpha = SoilPar.alpha;      % m^-1 Air entry value
    n = SoilPar.n;              % van Genuchten parameter
    m = 1 - 1 ./ SoilPar.n;     % van Genuchten parameter
    % Ksat = SoilPar.Ksat;        % m/s saturated hydraulic conductivity
    thetaS = SoilPar.thetaS;    % Saturated water content
    thetaR = SoilPar.thetaR;    % Residual water content

    rhow = 1000;
    g = 9.81;                   % m/s^2 accleration due to gravity
    Sw = 4e-10 .* rhow .* g;    % compressibility of water

    sEff = (1 + (alpha .* abs(h)) .^ n) .^ - m .* (h < 0) + (h >= 0);
    s = sEff + (thetaR ./ thetaS);
    c = m .* n .* alpha .* (alpha .* abs(h)) .^ (n - 1) .* ...
        (1 + (alpha .* abs(h)) .^ n) .^ (-m - 1) .* (h < 0) + 0 .* (h >= 0);
    c = c + Sw .* s;
    cMatr = sparse(diag(c, 0));
end

