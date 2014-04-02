function tPdf = TravelTimePdf(t, tFast, tSlow, sigma, mobileFrac)

    tPdf = mobileFrac * lognpdf(t, GetLognormalMu(tFast, sigma), sigma) + ...
        (1 - mobileFrac) * lognpdf(t, GetLognormalMu(tSlow, sigma), sigma);
    
    return
    
    function muX = GetLognormalMu(tMode, sigmaX)
        muX = log(tMode) + sigmaX * sigmaX;
    end
    
end