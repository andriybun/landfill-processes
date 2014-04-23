function FourierEstimateRainfall

    close all;

    CSV_FILENAME = ...
        '/home/andriybun/Workspace/tmp/Fellner/data.csv';

    num = dlmread(CSV_FILENAME, ',');
    num(isnan(num)) = 0;
    num = num(781:end, :);
    
    nT = size(num, 1);
    dt = 1;
    t = 1:dt:nT;
    
    qOut = num(:, 1);
    
    mu = -2;
    sigma = 1.5;
    
    %% Block for picking mu and sigma
    iStart = 230;
    nShow = 10;
    iShow = iStart:iStart+nShow;
    
    qResp = lognpdf(t(iShow) - iStart, mu, sigma) * dt;
    
    figure(1);
    plot(t(iShow), qOut(iShow));
    hold on;
    plot(t(iShow), 1e+4 * qResp, 'r');
    hold off;
    %% End block

    qRespF = fft(lognpdf(t', mu, sigma) * dt);
    qOutF = fft(qOut);
    rainDataF = qOutF ./ qRespF;
    rainDataEst = ifft(rainDataF);
    
    figure(2);
    plot(t, rainDataEst);
    hold on
    plot(t, qOut, 'r');
    hold off
    
end