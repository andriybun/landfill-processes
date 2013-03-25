function qTop = qBoundary(t)
    qTop = -0.02;
    tCond = 1;
    rainLength = 10;
    tCond = mod(t, 100) <= rainLength; % & ((t < 200) | (t >= 300));
    qTop = tCond * qTop;
end