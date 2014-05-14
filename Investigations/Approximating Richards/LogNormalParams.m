function out = LogNormalParams(theta)
    % 

    header = {'zref', 'theta', 'mu', 'sigma'};
    data = [
        -1.3, 0.117722, 9.895049, 0.630288;
        -1.2, 0.124992, 9.837246, 0.623807;
        -1.1, 0.135930, 9.762104, 0.613337;
        -1.0, 0.157630, 9.654392, 0.599806;
        -0.9, 0.184126, 9.501214, 0.588715;
        -0.8, 0.213004, 9.332375, 0.580947;
        -0.7, 0.241422, 9.153540, 0.575228;
        -0.6, 0.269278, 8.963267, 0.575914;
        -0.5, 0.296424, 8.760224, 0.588319;
        -0.4, 0.322634, 8.542915, 0.621239;
        -0.3, 0.347522, 8.312056, 0.672297;
        -0.2, 0.370333, 8.062729, 0.729063;
        -0.1, 0.389034, 7.675488, 0.770306;
        -0.05,0.394715, 7.085491, 0.795346
        ];
    
    mu = interp1(data(:, 2), data(:, 3), theta);
    sigma = interp1(data(:, 2), data(:, 4), theta);
    out = {mu, sigma};
    
%     fH = figure(1);
%     set(fH, 'Position', [500, 400, 700, 230]);
% %     title('Parameters of log-normal distribution as function of initial moisture content', ...
% %         'HorizontalAlignment', 'center');
%     h = subplot(1, 2, 1);
% %     set(h, 'position', [0, 0, 0.5, 0.8]);
%     plot(1 - data(:, 2) / 0.4, data(:, 3));
%     xlabel('Remaining storage capacity, %');
%     ylabel('mu');
%     xlim([0, 0.8]);
%     h = title('a)');
%     set(h, 'VerticalAlignment', 'bottom');
%     h = subplot(1, 2, 2);
% %     set(h, 'position', [0.5, 0, 0.5, 0.8]);
%     plot(1 - data(:, 2) / 0.4, data(:, 4));
%     xlabel('Remaining storage capacity, %');
%     ylabel('sigma');
%     xlim([0, 0.8]);
%     h = title('b)');
% %     set(h, 'VerticalAlignment', 'bottom');
% %     set(gcf, 'NextPlot', 'add');
% %     axes;
% %     strTitle = 'Parameters of log-normal distribution as function of initial moisture content';
% %     h = title(strTitle);
% %     set(gca, 'Visible', 'off');
% %     set(h, 'Visible', 'on');

end