function [qIn, qOut] = ComputeBoundaryFlux(qz, qx, ModelDim)
    % Top
    [qTopIn, qTopOut] = OneSideFlux(qz, 1, 1, ModelDim);
    % Bottom
    [qBottomIn, qBottomOut] = OneSideFlux(qz, 1, ModelDim.znin, ModelDim);
    % Left
    [qLeftIn, qLeftOut] = OneSideFlux(qx, 2, 1, ModelDim);
    % Right
    [qRightIn, qRightOut] = OneSideFlux(qx, 2, ModelDim.xnin, ModelDim);
    
    % Summarize
    qIn = qTopIn + qBottomIn + qLeftIn + qRightIn;
    qOut = qTopOut + qBottomOut + qLeftOut + qRightOut;
    
    return
    
    function [qInX, qOutX] = OneSideFlux(qX, dim, iInode, ModelDim)
        % dim: 1 - vertical, 2 - horizontal flux
        switch dim
            case 1
                maxIn = ModelDim.znin;
                dInSign = sign(ModelDim.dzin(1));
                qSide = qX(iInode, :);
            case 2
                maxIn = ModelDim.xnin;
                dInSign = sign(ModelDim.dxin(1));
                qSide = qX(:, iInode);
        end
        switch iInode
            case 1
                % Define alignment parameter to account for upstream/downstream boundaries
                boundSign = 1;
            case maxIn
                boundSign = -1;
        end

        % Flux aligned with axis (positive flux means influx, negative - outflux)
        qAlign = boundSign * dInSign * qSide;
        isInflux = (qAlign > 0);
        qInX = sum(qAlign(isInflux));
        qOutX = -sum(qAlign(~isInflux));
    end
end