function ModelDim = InitializeNodes(varName, varTop, varBottom, dVar, orientation)
%
% Function generates structure containing information about geometrical discretization of 
% investigated system (nodes / internodes).
% Function generates a list of field names corresponding to a variable characterizing one
% spatial dimension.
% Input parameters:
%   - varName           - is a name of axis for which discretization is generated
%   - varTop, varBottom - spatial boundaries of the system
%   - dVar              - spatial discretization along generated axis
%   - orientation       - one of the following: 'vertical' (default), 'horizontal'
% Fields generated are:
%   - positions of internodes (from 0 to varBottom with a step dVar)
%   - positions of nodes
%   - intervals between internodes
%   - intervals between nodes
%   - number of internodes
%   - number of nodes
%

    if nargin < 5
        orientation = 'vertical';
    end
    orientationId = find(strcmp({'vertical', 'horizontal'}, orientation), 1);

    % Templates for field names. Field names are generated using passed
    % varName.
    fnTemplates = {'%sin', '%sn', 'd%sin', 'd%sn', '%snin', '%snn',};
    vin  = sprintf(fnTemplates{1}, varName);
    vn   = sprintf(fnTemplates{2}, varName);
    dvin = sprintf(fnTemplates{3}, varName);
    dvn  = sprintf(fnTemplates{4}, varName);
    vnin = sprintf(fnTemplates{5}, varName);
    vnn  = sprintf(fnTemplates{6}, varName);
    
    % Initialize output structure
    ModelDim = struct();
    
    % Initialize its fields
    % Positions of internodes (equally spaced on interval [0, varBottom])
    switch (orientationId)
        case 1
            ModelDim.(vin) = (varTop:dVar:varBottom)';
        case 2
            ModelDim.(vin) = (varTop:dVar:varBottom);
        otherwise
            error('Initialization:SpatialDiscretization', 'Wrong orientation keyword');
    end
    % Number of internodes
    ModelDim.(vnin) = length(ModelDim.(vin));
    % Number of nodes
    ModelDim.(vnn) = ModelDim.(vnin) - 1;
    % Positions of nodes (in between internodes, first and last nodes are
    % aligned with first and last internodes respectively)
    ModelDim.(vn) = (ModelDim.(vin)(1:ModelDim.(vnin)-1) + ModelDim.(vin)(2:ModelDim.(vnin))) ./ 2;
    % Intervals between nodes
    ModelDim.(dvn) = ModelDim.(vn)(2:end) - ModelDim.(vn)(1:end-1);
    % Intervals between internodes
    ModelDim.(dvin) = ModelDim.(vin)(2:end) - ModelDim.(vin)(1:end-1);
end