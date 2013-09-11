function qTop = qBoundaryRainData(t)
    global rainData TimeParams
    iT = find(t >= TimeParams.t, 1, 'last');
    qTop = -rainData(iT);

%     qTop = -0.02;
%     tCond = 1;
%     rainLength = 10;
%     tCond = mod(t, 100) <= rainLength; % & ((t < 200) | (t >= 300));
%     qTop = tCond * qTop;
end