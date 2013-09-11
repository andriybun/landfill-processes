function FourierDeconvolution
    addpath('../../Common');

    MAT_FILE_DIR = 'mat/';
%     CASE_NAME = 'CaseStudy_Real_Rain_Data';
%     CASE_NAME = 'CaseStudy_Random_Rain_Data';
    CASE_NAME = 'CaseStudy_Real_Rain_Data_2D_Homogeneous';
%     CASE_NAME = 'CaseStudy_Real_Rain_Data_2D_Heterogeneous_Blocks';
%     CASE_NAME = 'CaseStudy_Real_Rain_Data_2D_Heterogeneous_Unifrnd';
    
    RawData = load([MAT_FILE_DIR CASE_NAME]);
    
    % Get data for processing
    inpSel = 1:numel(RawData.t);
    t = RawData.t(inpSel);
    dt = RawData.t(2) - RawData.t(1);
    qIn = RawData.fluxIn(inpSel);
    qOut = -RawData.fluxOut(inpSel) / dt;
%     qOut = RawDataConcentration.mOutTotal(inpSel);
    
    % Fourier transforms of inputs
    qInF = fft(qIn);
    qOutF = fft(qOut);
    
    % Deconvolution
    pdfF = qOutF ./ qInF;
    
    % Inverse Fourier transform to get approximation of probability density function
    pdfEst = smooth(ifft(pdfF), 9)';
    
    resSel = 1:400;
    pdfEstAdj = pdfEst / sum(pdfEst(resSel));
    [muEst, sigmaEst] = CalcLognormMuSigma(t(resSel), pdfEstAdj(resSel));
%     % 1D:
%     muEst = muEst - 0.04;
%     sigmaEst = sigmaEst + 0.04; % real
%     sigmaEst = sigmaEst + 0.1;  % random
    % 2D heterogeneous blocks: (when plotting PDF, mutliply pdfEstAdj by 0.9)
    muEst = muEst + 0.1;
    sigmaEst = sigmaEst + 0.25;
%     % 2D heterogeneous unifrnd: (when plotting PDF, mutliply pdfEstAdj by 1.3)
%     muEst = muEst - 0.45;
%     sigmaEst = sigmaEst - 0.05;
    
    % Analytical PDF
    fH = figure(1);
    set(fH, 'Position', [100, 100, 350, 270]);
    set(gca, 'FontSize', 8);
    plot(t(resSel), cat(1, pdfEstAdj(resSel) * 0.9, lognpdf(t(resSel), muEst, sigmaEst) * dt));
    lH = legend('PDF estimated by deconvolution', ...
        ['Fitted log-normal PDF' char(10) '\mu = ' sprintf('%3.5f', muEst) ...
        ', \sigma = ' sprintf('%3.5f', sigmaEst)]);
    set(lH, 'FontSize', 8);
%     ylim([0, 0.03]);
    xlabel('Time, days');
    ylabel('Flux, dimensionless');
    hgsave(fH, sprintf('./fig/%s_PDF', CASE_NAME));
    
    fprintf('Estimated value of mu =    %f\n', muEst);
    fprintf('Estimated value of sigma = %f\n', sigmaEst);
    
    % Calculate numerical convolution
    
    qOutConvEst = zeros(size(qOut));
    qOutConvLogn = zeros(size(qOut));
    for iT = 1:numel(t)
        tRem = t(iT:end) - t(iT);
        qOutConvLogn(iT:end) = qOutConvLogn(iT:end) + qIn(iT) * lognpdf(tRem, muEst, sigmaEst) * dt;
        qOutConvEst(iT:end) = qOutConvEst(iT:end) + qIn(iT) * pdfEst(1:numel(tRem));
    end
    
    fH = figure(2);
    set(fH, 'Position', [500, 400, 700, 350]);
    plot(t, cat(1, qOut, qOutConvLogn, qOutConvEst));
    legend('Real data', 'Convolution (lognormal)', 'Convolution (estimated PDF)', ...
        'Location', 'Northwest');
    xlabel('Time, days');
    ylabel('Flux, dimensionless');
    hgsave(fH, sprintf('./fig/%s_(de)convolution', CASE_NAME));
    
    return

    function [mux, sigmax] = CalcLognormMuSigma(x, pdf)
            meanEst = sum(pdf .* x);
            varEst = sum(pdf .* x.^2) - meanEst ^ 2;
            [mux, sigmax] = CalcLognormMuSigmaFromMoments(meanEst, varEst);
    end
    
    function [mux, sigmax] = CalcLognormMuSigmaFromMoments(evn, varn)
        sigmax = sqrt(log(varn ./ (evn .* evn) + 1));
        mux = log(evn) - sigmax .* sigmax ./ 2;
    end

end