function FourierConvolution
    addpath('../../Common/');
    addpath('../../Integrated Model/Data/');
    
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
    TimeParams = PrecipitationData.TimeParams;
    
    % Time parameters
    tEnd = TimeParams.t(end);   % 100;
    dt = TimeParams.dt;         % 1;
    t = TimeParams.t;           % 0:dt:tEnd;

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
    
    lognpdfVec = lognpdf(t, mu, sigma) * dt;

    tic
    
    rainDataF = fft(rainData);
    lognpdfVecF = fft(lognpdfVec);
    
    qOutF = rainDataF .* lognpdfVecF;
    qOut = ifft(qOutF);
    
    toc
    
    plot(t, qOut);
    
end