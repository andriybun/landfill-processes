function q = NettoFluxWater(t, h, SoilPar, ModelDim, BoundaryPar)
    % Function for calculating rates for a liquid flow in a spatial domain
    % qH are the water fluxes at all internodal points in the spatial domain;
    % t is time
    % T are the temperatures in all cells in the spatial domain
    % SoilPar contains the Model Parameters (thermal diffusivity of all
    % InterNodal points)
    % ModelDim contains the ModelDiscretization (coordinates of nodes and
    % internodes, deltas of nodal points and
    % internodal points)
    
    znn = ModelDim.znn;
    znin = ModelDim.znin;
    dzn = ModelDim.dzn;
    dzin = ModelDim.dzin;

    % Calculate Derived States (total head etc.)
    hIn(1, 1) = h(1, 1);
    hIn(2:znn, 1) = (dzin(1:znn-1, 1) .* h(1:znn-1, 1) + dzin(2:znn, 1) .* h(2:znn, 1)) ./ ...
        (dzin(1:znn-1, 1) + dzin(2:znn));
    hIn(znn+1, 1) = h(znn, 1);
    k = ComputeHydraulicConductivity(hIn, SoilPar, ModelDim);

    % Upper boundary zero flux
    q(1, 1) = BoundaryPar.qTop(t);
    q(2:znin-1, 1) = -k(2:znin-1, 1) .* ((h(2:znn, 1) - h(1:znn-1, 1)) ./ dzn(1:znn-1, 1) + 1);

    % Lower/Right boundary conditions
    % unit gradient (gravity drainage (zero pressure gradient) ...
    q(znin, 1) = BoundaryPar.kSurf .* (BoundaryPar.hAmb - h(znn, 1));
end