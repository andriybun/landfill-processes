function cOut = SoluteTransportAnalytic(x, t, u, d, ModelDim)
    nx = numel(x);
    nt = numel(t);

    checkInputs(x, u);
    checkInputs(x, d);
    
    u = -u;
    
    if ~(nx == 1 || nt == 1)
        x = repmat(x, [1, nt]);
        t = repmat(t, [nx, 1]);
        if numel(u) ~= 1
            u = repmat(u, [1, nt]);
        end
        if numel(d) ~= 1
            d = repmat(d, [1, nt]);
        end
    end
    
    x = -x;
    l = 1; % -ModelDim.zin(end);

    t0 = 500;
    Ci = 1;
    C0 = 0;
    isSolA = t < t0;
    
    solCase = 3;
    
     switch solCase
        %% From Weerts book:
        case 1
            % 4.1
            %  IC & BC:
            %   c(x, 0) = C_i;
            %   c(0, t) = C_0 (0 < t <= t_0) & 0 (t > t_0)
            %   dc/dx(inf, t) = 0
            cOut = zeros(size(t));
            cOut(isSolA) = Ci + (C0 - Ci) .* A1(x(isSolA), t(isSolA), l, u, d);
            cOut(~isSolA) = Ci + (C0 - Ci) .* A1(x(~isSolA), t(~isSolA), l, u, d) - ...
                C0 .* A1(x(~isSolA), t(~isSolA) - t0, l, u, d);
        case 2
            % 4.3
            %  IC & BC:
            %   c(x, 0) = C_i;
            %   c(0, t) = C_0 (0 < t <= t_0) & 0 (t > t_0)
            %   dc/dx(L, t) = 0
            cOut = zeros(size(t));
            cOut(isSolA) = Ci + (C0 - Ci) .* A3(x(isSolA), t(isSolA), l, u, d);
            cOut(~isSolA) = Ci + (C0 - Ci) .* A3(x(~isSolA), t(~isSolA), l, u, d) - ...
                C0 .* A3(x(~isSolA), t(~isSolA) - t0, l, u, d);

        %% Other solutions
        case 3
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
        case 4
            % IC: C = 1, 0 <= x < inf
            cOut = 0.5 - 0.5 * erf(0.5 * (-x + u .* t) ./ (sqrt(t .* d)));
            cOut = max(cOut, 0);
         case 5
             % IC: C = 0, C(0, t) = dirac(t)
             cOut = 0.1 * x .* exp(-(x - u .* t) .^ 2 ./ (4 .* d .* t)) ./ (2 .* sqrt(pi .* d .* t .^ 3));
    end
    
    return
    
    function checkInputs(spatial, var)
        if (numel(var) ~= 1) && (~isequal(size(var), size(spatial)))
            error('Wrong dimensions of input arguments');
        end
    
    end

    function res = A1(x, t, l, u, d)
        R = 1;
        res = 1 / 2 * erfc((R .* x - u .* t) ./ (2 .* sqrt(d .* R .* t))) + ...
            1 / 2 * exp(u .* x ./ d) .* erfc((R .* x + u .* t) ./ (2 .* sqrt(d .* R .* t)));
    end

    function res = A3(x, t, l, u, d)
        R = 1;
        res = 1 / 2 .* erfc((R .* x - u .* t) ./ (2 .* sqrt(d .* R .* t))) + ...
            1 / 2 .* exp(u .* x ./ d) .* erfc((R .* x + u .* t) ./ (2 .* sqrt(d .* R .* t))) + ...
            1 / 2 .* (2 + u .* (2 * l - x) ./ d + u .* u .* t ./ (d .* R)) * exp(u .* l ./ d) .* ...
            erfc((R .* (2 * l - x) + u .* t) ./ (2 .* sqrt(d .* R .* t))) - ...
            sqrt(u .* u .* t / (pi .* d .* R)) .* ...
            exp(u .* l ./ d - R ./ (4 .* d .* t) .* (2 * l - x + u .* t ./ R) .^ 2);
        
    end
end