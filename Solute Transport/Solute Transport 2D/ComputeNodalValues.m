function cNodes = ComputeNodalValues(zMark, xMark, cMark, ModelDim, StatType, ...
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

    % Nodal concentrations
    znNodes = ModelDim.znn;
    xnNodes = ModelDim.xnn;
    nSolutes = size(cMark, 2);
    cNodes = zeros(znNodes, xnNodes, nSolutes);
    
    for iZ = 1:znNodes
        for iX = 1:xnNodes
            isInNode = (zNodeIdx == iZ) & (xNodeIdx == iX);
            if any(isInNode)
                % Concentration at the node
                cNodes(iZ, iX, :) = StatType(cMark(isInNode, :));
            else
                % Zero concentration
                cNodes(iZ, iX, :) = 0;
            end
        end
    end

end
