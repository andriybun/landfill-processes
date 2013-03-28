function [mIn, mOut] = ComputeBoundaryMassFlux(qz, qx, cn, ModelDim, BoundaryPar, nSolutes)
    % Top
    [qTopIn, qTopOut] = OneSideFlux(qz, cn, 1, 1, ModelDim, BoundaryPar, nSolutes);
    % Bottom
    [qBottomIn, qBottomOut] = OneSideFlux(qz, cn, 1, ModelDim.znin, ModelDim, BoundaryPar, nSolutes);
    % Left
    [qLeftIn, qLeftOut] = OneSideFlux(qx, cn, 2, 1, ModelDim, BoundaryPar, nSolutes);
    % Right
    [qRightIn, qRightOut] = OneSideFlux(qx, cn, 2, ModelDim.xnin, ModelDim, BoundaryPar, nSolutes);
    
    % Summarize
    mIn = qTopIn + qBottomIn + qLeftIn + qRightIn;
    mOut = qTopOut + qBottomOut + qLeftOut + qRightOut;
    
    return
    
    function [mInX, mOutX] = OneSideFlux(qX, cX, dim, iInode, ModelDim, BoundaryPar, nSolutes)
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
                cBound = cX(min(iInode, ModelDim.znn), :, :);
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
                cBound = cX(:, min(iInode, ModelDim.xnn), :);
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
                mInX = sum(repmat(qAlign(isInflux), [1, 1, nSolutes]) .* ...
                    cBound(1.5 - 0.5 * boundSign, isInflux, :));
                mOutX = -sum(repmat(qAlign(~isInflux), [1, 1, nSolutes]) .* ...
                    cBound(1.5 + 0.5 * boundSign, ~isInflux, :));
            case 2
                mInX = sum(repmat(qAlign(isInflux), [1, 1, nSolutes]) .* ...
                    cBound(isInflux, 1.5 - 0.5 * boundSign, :));
                mOutX = -sum(repmat(qAlign(~isInflux), [1, 1, nSolutes]) .* ...
                    cBound(~isInflux, 1.5 + 0.5 * boundSign, :));
        end
        mInX = squeeze(mInX);
        mOutX = squeeze(mOutX);
    end
end