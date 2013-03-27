function qBnd = qBnd(t,z,x,BLoc,ModelDim)

%function will only be called from one boundary
qBnd = [];
nz = length(z);
nx = length(x);
switch lower(BLoc)
    case 'bottom',
        %top or bottom boundary
        qBnd(1,1:nx) = 0;
    case 'top'
        qBnd(1,1:nx) = -0.001;
    case 'left'
        qBnd(1:nz,1) = 0;
    case 'right'
        qBnd(1:nz,1) = 0;
end


end
