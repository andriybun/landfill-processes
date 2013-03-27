function [cNodes, sIdx] = ComputeNodalValues(zMark, xMark, cMark, ModelDim, StatType, ...
                                             zNodeIdx, xNodeIdx)
% Function:
%   Compute nodal values
%
% Purpose:
%   Function to calculate nodal values given values at a set of markers.
%   It uses markers with coordinates within half grid step from a given
%   node (between nearest internodes) to calculate values at any node.
%
% Parameters:
%   zMark, xMark - vectors of coordinates of markers.
%   cMark - vector of concentrations at markers. Must be of the same
%           size as z.
%   ModelDim - structure with dimensions of model.
%   StatType - handle of statistical function (e.g. sum, mean, max etc.)
%   (optional) zNodeIdx, xNodeIdx - indices of nodes to which markers belong
%
% Return:
%   cNodes - concentration of values at given nodes.
%   sIdx   - indices of elements in zMark while sorted.

    if nargin < 4
        StatType = @mean;
    end

    % Sort markers
    [zMarkSorted, sIdx] = sort(-zMark);
    xMarkSorted = xMark(sIdx, :);
    
    %% Alternative using average
    cMarkSorted = cMark(sIdx, :);

    % Nodal concentrations
    znNodes = ModelDim.znn;
    xnNodes = ModelDim.xnn;
    nSolutes = size(cMark, 2);
    cNodes = zeros(znNodes, xnNodes, nSolutes);
    
    % Index of nearest node to marker
    if nargin < 6
        zNodeIdx = interp1(ModelDim.zn, 1:ModelDim.znn, -zMarkSorted, 'nearest', 'extrap');
        xNodeIdx = interp1(ModelDim.xn, 1:ModelDim.xnn, -xMarkSorted, 'nearest', 'extrap');
    end
    
    for iZ = 1:znNodes
        for iX = 1:xnNodes
            isInNode = (zNodeIdx == iZ) & (xNodeIdx == iX);
            if any(isInNode)
                % Concentration at the node
                cNodes(iZ, iX, :) = StatType(cMarkSorted(isInNode, :));
            else
                % Zero concentration
                cNodes(iZ, iX, :) = 0;
            end
        end
    end

end
