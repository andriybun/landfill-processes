function nf = NettoFluxConcNodes(c, d, theta, ModelDim)
    % Input parameters:
    %   c       - current concentration of solute
    %   d       - diffusion coefficient

    qIn = SoluteFlow1d(c, d, theta, ModelDim);
    nNodes = ModelDim.znn;
    dzin = ModelDim.dzin;
    
    nSolutes = size(c, 2);
    
    nf = zeros(nNodes, nSolutes);
    allNodes = 1:nNodes;
    
    for soluteIdx = 1:nSolutes
        nf(allNodes, soluteIdx) = ...
            -(qIn(allNodes + 1, soluteIdx) - qIn(allNodes, soluteIdx)) ./ (dzin(allNodes));
    end
    
    nf = nf ./ repmat(theta, [1, nSolutes]);
    
    return
    
    function qIn = SoluteFlow1d(c, d, theta, ModelDim)
        nInodes = ModelDim.znin;
        dzn = ModelDim.dzn;
        
        nSolutes = size(c, 2);
        
        % Diffusive Flux
        dIn = DistributeValue(d, ModelDim.zin);
        thetaIn = zeros(nInodes, 1);
        
        % Internodal flux
        qIn = zeros(nInodes, 1);
        internInodes = 2:nInodes - 1;
        
        thetaIn(internInodes) = (theta(internInodes - 1) + theta(internInodes)) / 2;
        
        for xSoluteIdx = 1:nSolutes
            qIn(internInodes, xSoluteIdx) = ...
                -dIn(internInodes, xSoluteIdx) .* thetaIn(internInodes) .* ...
                (c(internInodes, xSoluteIdx) - c(internInodes - 1, xSoluteIdx)) ./ ...
                dzn(internInodes - 1);
        end        
        
    end

    function valDistr = DistributeValue(val, nodes)
        valDistr =  ones(size(nodes)) * val;
    end
end