function netInf = NetInfiltration(rf, ev, alpha)
    % Function to calclulate net infiltration from rainfall and evapotranspiration
    
    %% Option 1: limit net infiltration to zero + evaporate water from previous time step
%     netInf = rf - ev;
%     netInf(1) = max(0, netInf(1));
%     isNegative = (netInf < 0);
%     isNegativeShift = cat(1, isNegative(2:end), false);
%     netInf(isNegativeShift) = netInf(isNegativeShift) + netInf(isNegative);
%     netInf(isNegative) = 0;
%     netInf(isNegativeShift) = max(0, netInf(isNegativeShift));
    
    %% Option 2: subtract
    if nargin < 3
        alpha = 1;
    end
    netInf = rf - alpha * ev;
end