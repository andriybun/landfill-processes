function SoilPar = InitializeSoilProperties(kSat, ModelDim, mobileFraction)
    hCap = 0.5;
    SoilPar.alpha = 1 / hCap;
    SoilPar.n = 1.9;
    
    SoilPar.thetaR = 0.04;
    SoilPar.thetaS = 0.4;
    
    if nargin >= 3
        SoilPar.thetaR = SoilPar.thetaR * mobileFraction;
        SoilPar.thetaS = SoilPar.thetaS * mobileFraction;
    end
    
    SoilPar.kSat = kSat;
    
    SoilPar.theta = 0.3;
end