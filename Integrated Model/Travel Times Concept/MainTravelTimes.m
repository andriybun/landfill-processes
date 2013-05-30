function MainTravelTimes
    addpath('../../Common/');
    addpath('../Data/');
    
    % Tolerance
    EPSILON = 1e-10;
    NUM_SIGMAS = 6;
    
    % Dimensions
    zTop = 0;
    zBottom = -1;
    dz = 0.05;
    ModelDim = InitializeNodes('z', zBottom, zTop, dz);
    
    % Load: rainData, TimeParams, StartDate
    PrecipitationData = load('precipitation');
    % Precipitation in meters per time interval dt
    rainData = 4 * PrecipitationData.rainData;
%     rainData = 1e-2 * ones(size(rainData));
    TimeParams = PrecipitationData.TimeParams;
    
    % Time parameters
    tEnd = TimeParams.t(end);   % 100;
    dt = TimeParams.dt;         % 1;
    t = TimeParams.t;           % 0:dt:tEnd;
    nT = TimeParams.numIntervals;

    % Log-normal parameters
    mu = 0;
    sigma = 0.6;

    % Other geometrical
    zLength = abs(ModelDim.zin(1) - ModelDim.zin(end));
    % Pore volume
    theta = 0.4;
    pv = zLength * theta;
    
    % Decay rate
    lambda = 1e-2;
    % Exchange rate between mobile-immobile phases
    kExch = 1e-0;
    
    % Initial concentration (immobile, mobile phases)
    cIni = [1; 0];
    % Initial mass of solute
    mIni = cIni * pv;
    
%     TimeParams.maxDays = 30;
    nT = TimeParams.maxDays * TimeParams.intervalsPerDay;
    t = t(1:nT);
    
    qOutTotal = zeros(1, nT);
    mOutTotal = zeros(1, nT);
    cRemaining = zeros(2, nT + 1);
    cRemaining(:, 1) = cIni;
    
    profile on
    tic
    
    % mRemaining = mIni;
    for iT = 1:nT
        tOffset = t(iT);
        tAfter = t(iT:nT) - tOffset;
        
        % Calculate boundaries of log-normally distributed travel times (to save computational
        % resources we calculate only for those times, when something comes out of the drainage
        % system but not until the end of the period).
        tBouds = LogNormalBounds(mu, sigma, NUM_SIGMAS);
        % Select only those time steps that are affected by current injection
        iCalcT = (tAfter <= tBouds(2));
        nCalcT = sum(iCalcT);
        iTend = iT + nCalcT - 1;
        tAfter = tAfter(iCalcT);
        
        % Every input impulse of water will cause (log-normal) response at the outlet. All the
        % outflow during a given time step is considered as a particle with unique travel time.
        qOutAfter = rainData(iT) * lognpdf(tAfter, mu, sigma) * dt;
        % We integrate volumes of all the particles flowing out at the same time intervals to 
        % obtain Leachate volume flux.
        qOutTotal(iT:iTend) = qOutTotal(iT:iTend) + qOutAfter;
        % The longer particle resides inside the landfill, the bigger oncentrations of solutes
        % will be.
        cOutAfter = OutConcentration(tAfter, cRemaining(:, iT), kExch, lambda);
        
        % Calculate the remaining mass of solute
        mRemaining = cOutAfter(:, 1) * pv;
        % We integrate masses of solutes in all particles.
        mOutTotal(iT:iTend) = mOutTotal(iT:iTend) + qOutAfter .* cOutAfter(2, :);
        mRemaining(2) = mRemaining(2) - mOutTotal(iT);
        
        % Calculate the remaining concentrations
        cRemaining(:, iT + 1) = mRemaining / pv;
        
        if mod(iT, 87) == 0
            fprintf('%5.3f%% complete\n', iT / nT * 100);
        end
    end

    toc
    profile off
    profile viewer
    
    %% Error check
    mEnd = cRemaining(:, end) * pv;
    if ~RealEq(sum(mIni) - sum(mOutTotal), sum(mEnd), EPSILON)
        warning('ResultCheck:MassBalanceError', 'Absolute error is too high: err = %3.2e', ...
            abs(abs(sum(mIni) - sum(mOutTotal) - sum(mEnd))));
    end
    %% End error check
    
%     return
    
    
    
    %% Plotting
    figH = figure(1);
    [axH, lH1, lH2] = plotyy(t, rainData(1:nT), t, sum(cRemaining(:, 1:nT), 1) / 2);
    xlabel('time [days]');
    set(get(axH(1), 'ylabel'), 'string', 'precipitation [m/hour]');
    set(get(axH(2), 'ylabel'), 'string', 'emission potential [m^3/m^3]');
    hgsave(sprintf('fig/oug_flux_c_rem_%d_days_lambda_%4.3f.fig', ...
        ceil(nT / TimeParams.intervalsPerDay), lambda));
%     print(figH, '-dpng', '-r0', sprintf('fig/oug_flux_c_rem_%d_days_lambda_%4.3f.png', ...
%         ceil(nT / TimeParams.intervalsPerDay), lambda));

    figH = figure(2);
    figPos = [100, 100, 500, 250];
    set(figH, 'Position', figPos);
    plot(t, rainData(1:nT), 'b');
    hold on;
    plot(t, qOutTotal(1:nT), 'r');
    hold off;
    xlabel('time [days]');
    ylabel('flux [m/hour]')
    legend({'Precipitation', 'Leachate volume flux'}, 'Location', 'NorthEast');
    hgsave(sprintf('fig/precip_leachate_flux_%d_days_lambda_%4.3f.fig', ...
        ceil(nT / TimeParams.intervalsPerDay), lambda));
%     print(figH, '-dpng', '-r0', sprintf('fig/precip_leachate_flux_%d_days_lambda_%4.3f.png', ...
%         ceil(nT / TimeParams.intervalsPerDay), lambda));
    
    figH = figure(3);
    figPos = [200, 200, 500, 300];
    set(figH, 'Position', figPos);
    [axH, lH1, lH2] = plotyy(t, qOutTotal(1:nT), t, mOutTotal(1:nT) ./ qOutTotal(1:nT));
    legend({'Leachate volume flux', 'Leachate concentration'}, 'Location', 'NorthEast');
    xlabel('time [days]');
    set(get(axH(1), 'ylabel'), 'string', 'out flux [m/hour]');
    set(get(axH(2), 'ylabel'), 'string', 'out concentration [m^3/m^3]');
    hgsave(sprintf('fig/concentr_leachate_flux_%d_days_lambda_%4.3f.fig', ...
        ceil(nT / TimeParams.intervalsPerDay), lambda));
%     print(figH, '-dpng', '-r0', sprintf('fig/concentr_leachate_flux_%d_days_lambda_%4.3f.png', ...
%         ceil(nT / TimeParams.intervalsPerDay), lambda));
    
    figH = figure(4);
    figPos = [300, 300, 500, 300];
    set(figH, 'Position', figPos);
    plot(t, mOutTotal(1:nT) ./ qOutTotal(1:nT), 'r');
    hold on;
    plot(t, sum(cRemaining(:, 1:nT), 1) / 2);
    hold off;
    xlabel('time [days]');
    ylabel('concentration [m^3/m^3]')
    legend({'Leachate concentration', 'Emission potential'}, 'Location', 'SouthEast'); % 'SouthWest'
    hgsave(sprintf('fig/concentr_c_rem_%d_days_lambda_%4.3f.fig', ...
        ceil(nT / TimeParams.intervalsPerDay), lambda));
%     print(figH, '-dpng', '-r0', sprintf('fig/concentr_c_rem_%d_days_lambda_%4.3f.png', ...
%         ceil(nT / TimeParams.intervalsPerDay), lambda));
    
    return
    
    function cOut = OutConcentration(tau, cRemaining, kExch, lambda)
%         cOut = cRemaining * (1 - exp(-lambda * tau));
        if (tau(1) == tau(end))
            cOut = cRemaining;
        else
            tRange = [tau(1), tau(end)];
            [tauOde, cOde] = ode45(@(tauX, cX) Dc(tauX, cX, kExch, lambda), tRange, cRemaining);
            cOut = interp1(tauOde, cOde, tau)';
        end
        
        return
        
        function dcdt = Dc(xTau, xC, kExch, lambda)
            dcdt = nan(2, 1);
            dcdt(1) = (-kExch - lambda) * xC(1) + kExch * xC(2);
            dcdt(2) = kExch * xC(1) - kExch * xC(2);
        end
    end
end