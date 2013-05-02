function Main
    addpath('../../Common/');
    
    zTop = 0;
    zBottom = -1;
    dz = 0.05;
    ModelDim = InitializeNodes('z', zBottom, zTop, dz);
    
    tEnd = 100;
    dt = 1;
    t = 0:dt:tEnd;

%     seed = 1;
%     rng(seed);
    
    u = -5e-2;
    nT = numel(t);
    scale = rand(1, nT);
    
    % Log-normal parameters
    mu = -2;
    sigma = 0.6;
    
    qOut = zeros(1, nT);
    for iT = 1:nT
        tOffset = t(iT);
        qOut(iT:nT) = qOut(iT:nT) + scale(iT) * lognpdf(t(iT:nT) - tOffset, mu, sigma);
    end
    figure(2);
    plotyy(t, scale, t, qOut);
    
    uScaled = u * scale;
    
    l = 1;
    d = 1e-3;
    
    % Baseline concentration
    cBase = SoluteTransport(t, ModelDim.zn, u, l, d);
    cScaled = ScaleTime(@SoluteTransport, t, ModelDim.zn, uScaled, l, d);
    
    % Plot
    selectedNodes = 1:5:ModelDim.znn;
    doDisplayScaled = true;
    DisplayPlot(t, cBase, cScaled, ModelDim, d, u, selectedNodes, doDisplayScaled);
end