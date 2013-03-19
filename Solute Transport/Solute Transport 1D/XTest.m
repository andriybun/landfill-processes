function XTest
    addpath('PhysicalProcesses');

    d = 1e-3;
    
    zTop = 0;
    zBottom = -1;
    dz = -0.05;
    ModelDim = InitializeNodes('z', zTop, zBottom, dz);
    nNodes = ModelDim.znn;
    
    thetaNodesX = ones(nNodes, 1);
    cNodesSol = ones(nNodes, 1);
    cNodesSol(1:3) = [0.1, 0.3, 0.8];

    tRangeOde = [0, 500];
    
    optionsOde = odeset('RelTol', 1e-5, 'AbsTol', 1e-3);
    [~, cNodesDiff] = ode15s(...
        @(t, cN) NettoFluxConcNodes(cN, ...
        d, ...
        thetaNodesX, ...
        ModelDim), ...
        tRangeOde, ...
        cNodesSol, ...
        optionsOde);
    
    cNext = cNodesDiff(end, :)';

    mesh(cNodesDiff);
end