function k = ComputeHydraulicConductivity(h, SoilPar, ModelDim)
    znn = ModelDim.znn;
    znin = ModelDim.znin;
    dzin = ModelDim.dzin;
    
    alpha = SoilPar.alpha;
    n = SoilPar.n;
    m = 1 - 1 ./ n;
    kSat = SoilPar.kSat;
    
    % Calculate internodal pressures
    hIn(1, 1) = h(1, 1);
    hIn(2:znn, 1) = (dzin(1:znn-1, 1) .* h(1:znn-1, 1) + dzin(2:znn, 1) .* h(2:znn, 1)) ./ ...
        (dzin(1:znn-1, 1) + dzin(2:znn));
    hIn(znn+1, 1) = h(znn, 1);
    
    allInodes = 1:znin;
    sEffIn(allInodes, 1) = (1 + (alpha .* abs(hIn)) .^ n).^ - m .* (hIn < 0) + (hIn >= 0);
    krw = sEffIn .^ 3;
    k = kSat .* krw;
      
end

