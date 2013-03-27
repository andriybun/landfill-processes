function hBound = hBound(t,z,x,ModelDim)

%function will only be called from one boundary
nz = length(z);
nx = length(x);
if nz == 1,
    %top or bottom boundary
    if z == ModelDim.zIN(1),
        hBound(1:nx) = -3;
    end,
    if z == ModelDim.zIN(end),
        hBound (1:nx) = -0.5;
    end
end

if nx == 1
    %top o bottomboundary
    if x == ModelDim.xIN(1),
        hBound(1:nz) = -3;
    end
    if x == ModelDim.xIN(end),
        hBound (1:nz) = -0.5;
    end
end

end
