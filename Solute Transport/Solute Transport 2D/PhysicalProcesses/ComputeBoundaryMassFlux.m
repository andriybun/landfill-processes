function [qIn, qOut] = ComputeBoundaryMassFlux(qz, qx, cn, ModelDim, BoundaryPar)
    % Top
    [qTopIn, qTopOut] = OneSideFlux(qz, cn, 1, 1, ModelDim, BoundaryPar);
    % Bottom
    [qBottomIn, qBottomOut] = OneSideFlux(qz, cn, 1, ModelDim.znin, ModelDim, BoundaryPar);
    % Left
    [qLeftIn, qLeftOut] = OneSideFlux(qx, cn, 2, 1, ModelDim, BoundaryPar);
    % Right
    [qRightIn, qRightOut] = OneSideFlux(qx, cn, 2, ModelDim.xnin, ModelDim, BoundaryPar);
    
    % Summarize
    qIn = qTopIn + qBottomIn + qLeftIn + qRightIn;
    qOut = qTopOut + qBottomOut + qLeftOut + qRightOut;
    
    return
    
    function [qInX, qOutX] = OneSideFlux(qX, cX, dim, iInode, ModelDim, BoundaryPar)
        cFieldList = {'cTop', 'cBottom', 'cLeft', 'cRight'};
        
        % dim: 1 - vertical, 2 - horizontal flux
        switch dim
            case 1
                % Maximum available internode index
                maxIn = ModelDim.znin;
                % Determine what is the positive direction along chosen axis
                dInSign = sign(ModelDim.dzin(1));
                % Flux across boundary
                qBound = qX(iInode, :);
                % Concentration of nodes nearest to boundary
                cBound = cX(min(iInode, ModelDim.znn), :);
                % Define field name of boundary concentration to be used
                cFieldList = cFieldList([1, 2]);
            case 2
                % Maximum available internode index
                maxIn = ModelDim.xnin;
                % Determine what is the positive direction along chosen axis
                dInSign = sign(ModelDim.dxin(1));
                % Flux across boundary
                qBound = qX(:, iInode);
                % Concentration of nodes nearest to boundary
                cBound = cX(:, min(iInode, ModelDim.xnn));
                % Define field name of boundary concentration to be used
                cFieldList = cFieldList([3, 4]);
        end
        switch iInode
            case 1
                % Define alignment parameter to account for upstream/downstream boundaries
                boundSign = 1;
                % Define field name of boundary concentration to be used
                cField = cFieldList{1};
                % Concatenate concentrations of nodes at boundary with boundary concentrations
                cBound = cat(dim, BoundaryPar.(cField), cBound);
            case maxIn
                boundSign = -1;
                cField = cFieldList{2};
                cBound = cat(dim, cBound, BoundaryPar.(cField));
        end

        % Flux aligned with axis (positive flux means influx, negative - outflux)
        qAlign = boundSign * dInSign * qBound;
        isInflux = (qAlign > 0);
        switch dim
            case 1
                qInX = sum(qAlign(isInflux) .* cBound(1.5 - 0.5 * boundSign, isInflux));
                qOutX = -sum(qAlign(~isInflux) .* cBound(1.5 + 0.5 * boundSign, ~isInflux));
            case 2
                qInX = sum(qAlign(isInflux) .* cBound(isInflux, 1.5 - 0.5 * boundSign));
                qOutX = -sum(qAlign(~isInflux) .* cBound(~isInflux, 1.5 + 0.5 * boundSign));
        end
        
    end
end