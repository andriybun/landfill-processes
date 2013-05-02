function y = ScaleTime(f_handle, t, x, u, varargin)
    y = nan(numel(x), numel(t));
    udt = 0;
    u0 = u(1);

    % step t = 0
    y(:, 1) = f_handle(t(1), x, u0, varargin{:});

    % loop over t
    for idx = 2:numel(t)
        ratio = u(idx) / u0;
        udt = udt + (t(idx) - t(idx-1)) * ratio;
        y(:, idx) = f_handle(udt, x, u0, varargin{:});  % 1, numel(t)
    end
end