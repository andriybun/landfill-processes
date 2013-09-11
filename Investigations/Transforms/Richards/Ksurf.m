function [Kval,hAmb] = Ksurf(t,z,x,BLoc,ModelDim)
% Robbins condition at left and right domain edges (not implemented)
Kval = [];
hAmb = [];
nz = length(z);
nx = length(x);
switch BLoc,
    case 'bottom'
        Kval(1,1:nx) = 500;
        hAmb(1,1:nx) = -1.1;
    case 'top'
        Kval(1,1:nx) = 500;
        hAmb (1,1:nx) = 0;
    case 'left'
        Kval(1:nz,1) = 0.0;
        hAmb(1:nz,1) = -3;
    case 'right'
        Kval(1:nz,1) = 0;
        hAmb(1:nz,1) = -0.5;
end
