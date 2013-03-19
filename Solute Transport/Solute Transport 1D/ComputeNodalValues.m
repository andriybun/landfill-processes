function [cNodes, sIdx] = ComputeNodalValues(zMark, cMark, ModelDim, StatType, zNodeIdx)
% Function:
%   Compute nodal values
%
% Purpose:
%   Function to calculate nodal values given values at a set of markers.
%   It uses markers with coordinates within half grid step from a given
%   node (between nearest internodes) to calculate values at any node.
%
% Parameters:
%   zMark - vector of coordinates of markers.
%   cMark - vector of concentrations at markers. Must be of the same
%           size as z.
%   zNode - coordinates of nodes.
%
% Return:
%   cNodes - concentration of values at given nodes.
%   sIdx   - indices of elements in zMark while sorted.

    if nargin < 4
        StatType = @mean;
    end

    % Sort markers
    [zMarkSorted, sIdx] = sort(-zMark);

%     %% Alternative calculation using interp1:
%     cNodes = interp1(zMark(sIdx), cMark(sIdx), ModelDim.zn, INTERPOLATION_METHOD, 'extrap');
    
    %% Alternative using average
    cMarkSorted = cMark(sIdx, :);

    % Boundaries of cells to consider
    zInDown = -ModelDim.zin(1:end-1);
    zInUp = -ModelDim.zin(2:end);

    % Nodal concentrations
    nNodes = ModelDim.znn;
    nSolutes = size(cMark, 2);
    cNodes = zeros(nNodes, nSolutes);
    
    % Index of nearest node to marker
    if nargin < 5
        zNodeIdx = interp1(ModelDim.zn, 1:ModelDim.znn, -zMarkSorted, 'nearest', 'extrap');
    end
    
    for idx = 1:nNodes
        isInNode = (zNodeIdx == idx);
        if any(isInNode)
            % Concentration at the node
            cNodes(idx, :) = StatType(cMarkSorted(isInNode, :));
        else
            % Zero concentration
            cNodes(idx, :) = 0;
        end
    end

end
