function [States] = ComputeWaterStates(t, h, SoilPar, ModelDim, BoundaryPar)
    alpha = SoilPar.alpha;
    n = SoilPar.n;
    m = 1 - 1 ./ n;
    thetaR = SoilPar.thetaR;
    thetaS = SoilPar.thetaS;
    
    sEff = (1 + (alpha .* abs(h)).^n).^-m .* (h < 0) + (h >= 0);
    theta = thetaR + sEff .* (thetaS - thetaR);
    
    q = NettoFluxWater(t, h, SoilPar, ModelDim, BoundaryPar);
    
    States = struct();
    States.sEff = sEff;
    States.theta = theta;
    States.q = q;
end