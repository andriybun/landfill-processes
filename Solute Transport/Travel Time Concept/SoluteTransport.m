function cOut = SoluteTransport(t, x, u, l, d)
    nx = numel(x);
    nt = numel(t);

    checkInputs(t, u);
    checkInputs(x, d);
    
    u = -u;
    x = -x;
    
    if ~(nx == 1 || nt == 1)
        x = repmat(x, [1, nt]);
        t = repmat(t, [nx, 1]);
        if numel(u) ~= 1
            u = repmat(u, [nx, 1]);
        end
        if numel(d) ~= 1
            d = repmat(d, [1, nt]);
        end
    end
    
    % Analytical solution for the initial value problem:
    % IC: C = 1, 0 <= x <= l
    doReflect = true;
    if doReflect
        x = x + 1e-14;
        cOut = -(1 ./ 2) .* erf((1 ./ 2) .* (-x + u .* t) ./ (sqrt(t) .* sqrt(d))) + ...
            (1 ./ 2) .* erf((1 ./ 2) .* (-x + l + u .* t) ./ (sqrt(t) .* sqrt(d))) - ...
            (1 ./ 2) .* erf((1 ./ 2) .* (x + u .* t) ./ (sqrt(t) .* sqrt(d))) + ...
            (1 ./ 2) .* erf((1 ./ 2) .* (x + l + u .* t) ./ (sqrt(t) .* sqrt(d)));
    else
        cOut = -0.5 .* erf(0.5 .* (-x + u .* t) ./ sqrt(d .* t)) + ...
            0.5 .* erf(0.5 .* (-x + l + u .* t) ./ sqrt(d .* t));
    end

    return
    
    function checkInputs(spatial, var)
        if (numel(var) ~= 1) && (~isequal(size(var), size(spatial)))
            error('Wrong dimensions of input arguments');
        end
    end
end