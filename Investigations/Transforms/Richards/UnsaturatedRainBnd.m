function qBnd = UnsaturatedRainBnd(t,z,x,BLoc,ModelDim)

    global rainData TimeParams;

    %function will only be called from one boundary
    qBnd = [];
    nz = length(z);
    nx = length(x);
    switch lower(BLoc)
        case 'bottom',
            qBnd(1,1:nx) = 0;
        case 'top'
            iT = find(TimeParams.t < t + 1e-6, 1, 'last');
            qBnd(1,1:nx) = -rainData(iT) * TimeParams.dt; % (TimeParams.t(iT+1) - t);
        case 'left'
            qBnd(1:nz,1) = 0;
        case 'right'
            qBnd(1:nz,1) = 0;
    end


end