function y = ComputeY(t, h, SoilPar, ModelDim, BoundaryPar)
    
    znn = ModelDim.znn;
    allN = 1:znn;
    znin = ModelDim.znin;
    allIn = 1:znin;
    
    dzn = ModelDim.dzn;
    dzin = ModelDim.dzin;
    
    k = ComputeHydraulicConductivity(h, SoilPar, ModelDim);
    kSurf = BoundaryPar.kSurf;

    y(1, 1) = k(2, 1) / dzin(1) + BoundaryPar.qTop(t) ./ dzin(1);
    y(2:znn-1, 1) = (-k(2:znn-1, 1) + k(3:znn, 1)) ./ dzin(2:znn-1, 1);
    y(znn, 1) = -kSurf ./ dzin(znn) .* BoundaryPar.hAmb - k(znn) ./ dzin(znn);
end



