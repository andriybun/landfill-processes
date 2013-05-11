function MainTravelTimes
    addpath('../../Common/');
    addpath('../Data/');
    
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
    nT = TimeParams.num_intervals;

    % Log-normal parameters
    mu = 0;
    sigma = 0.6;

    % Other geometrical
    zLength = abs(ModelDim.zin(1) - ModelDim.zin(end));
    % Pore volume
    theta = 0.4;
    pv = zLength * theta;
    
    % Decay rate
    lambda = 1e-0;
    
    % Initial concentration
    cIni = 1;
    % Initial mass of solute
    mIni = pv * cIni;
    
    nT = 720;
    t = t(1:nT);
    
    qOutTotal = zeros(1, nT);
    mOutTotal = zeros(1, nT);
    cRemaining = zeros(1, nT + 1);
    cRemaining(1) = cIni;
    mRemaining = mIni;
    for iT = 1:nT
        tOffset = t(iT);
        tAfter = t(iT:nT) - tOffset;
        % Every input impulse of water will cause (log-normal) response at the outlet. All the
        % outflow during a given time step is considered as a particle with unique travel time.
        qOutAfter = rainData(iT) * lognpdf(tAfter, mu, sigma) * dt;
        % We integrate volumes of all the particles flowing out at the same time intervals to 
        % obtain leachate load.
        qOutTotal(iT:nT) = qOutTotal(iT:nT) + qOutAfter;
        % The longer particle resides inside the landfill, the bigger oncentrations of solutes
        % will be.
        cOutAfter = OutConcentration(tAfter, cRemaining(iT), lambda);
        % We integrate masses of solutes in all particles.
        mOutTotal(iT:nT) = mOutTotal(iT:nT) + qOutAfter .* cOutAfter;
        % Calculate the remaining mass of solute
        mRemaining = mRemaining - mOutTotal(iT);
        cRemaining(iT + 1) = mRemaining / pv;
    end

    figure(1);
    plotyy(t, rainData(1:nT), t, cRemaining(1:nT));

    figure(2);
    plot(t, rainData(1:nT), 'b');
    hold on;
    plot(t, qOutTotal(1:nT), 'r');
    hold off;
    
    figure(3);
    plotyy(t, qOutTotal(1:nT), t, mOutTotal(1:nT) ./ qOutTotal(1:nT));

    return
    
    function cOut = OutConcentration(t, cRemaining, lambda)
        cOut = cRemaining * (1 - exp(-lambda * t));
    end
end