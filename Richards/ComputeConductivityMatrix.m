function kMatrSparse = ComputeConductivityMatrix(h, SoilPar, ModelDim, BoundaryPar)
    
    znn = ModelDim.znn;
    allN = 1:znn;
    znin = ModelDim.znin;
    allIn = 1:znin;
    
    dzn = ModelDim.dzn;
    dzin = ModelDim.dzin;
    
    k = ComputeHydraulicConductivity(h, SoilPar, ModelDim);
    kSurf = BoundaryPar.kSurf;
     
    a(1, 1) = 0;
    b(1, 1) = -k(2, 1) ./ (dzin(1, 1) .* dzn(1, 1));
    c(1, 1) = k(2, 1) ./ (dzin(1, 1) .* dzn(1, 1));
    
    % Middle nodes
    a(2:znn-1, 1) = k(2:znn-1, 1) ./ (dzin(2:znn-1) .* dzn(1:znn-2));
    b(2:znn-1, 1) = -k(2:znn-1, 1) ./ (dzin(2:znn-1, 1) .* dzn(1:znn-2)) - ...
        k(3:znn, 1) ./ (dzin(2:znn-1, 1) .* dzn(2:znn-1, 1));
    c(2:znn-1, 1) = k(3:znn, 1) ./ (dzin(2:znn-1) .* dzn(2:znn-1, 1));
    
    % Lower boundary
    a(znn) = k(znn) ./ (dzin(znn) .* dzn(znn-1));
    b(znn) = -k(znn) ./ (dzin(znn) .* dzn(znn-1)) + kSurf ./ dzin(znn);
    c(znn) = 0;
    
    kMatr = diag(a(2:znn), -1) + diag(b, 0) + diag(c(1:znn-1), +1);
    kMatrSparse = sparse(kMatr);
end



