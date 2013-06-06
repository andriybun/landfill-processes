function ShowPlots(qOutTotal, mOutTotal, cRemaining, rainData, lambda, TimeParams)

    nT = TimeParams.maxDays * TimeParams.intervalsPerDay;
    t = TimeParams.t(1:nT);

    %% Prepare data structures to be plotted:
    TInfo = struct();
    TInfo.data = t;
    TInfo.axisLabel = 'time [days]';
    
    PrecipInfo = struct();
    PrecipInfo.data = rainData(1:nT);
    PrecipInfo.name = 'Precipitation';
    PrecipInfo.axisLabel = 'precipitation [m/hour]';
    PrecipInfo.color = [0.4, 0.7, 1];
    
    EmissionPotentialInfo = struct();
    EmissionPotentialInfo.data = sum(cRemaining(:, 1:nT), 1) / 2;
    EmissionPotentialInfo.name = 'Emission Potential';
    EmissionPotentialInfo.axisLabel = 'concentration [m^3/m^3]';
    EmissionPotentialInfo.color = [0.4, 0.4, 0.4];
    
    LeachateFluxInfo = struct();
    LeachateFluxInfo.data = qOutTotal(1:nT);
    LeachateFluxInfo.name = 'Leachate volume flux';
    LeachateFluxInfo.axisLabel = 'out flux [m/hour]';
    LeachateFluxInfo.color = [0.55, 0.25, 0.08];
    
    LeachateConcentrationInfo = struct();
    LeachateConcentrationInfo.data = mOutTotal(1:nT) ./ qOutTotal(1:nT);
    LeachateConcentrationInfo.name = 'Leachate concentration';
    LeachateConcentrationInfo.axisLabel = 'concentration [m^3/m^3]';
    LeachateConcentrationInfo.color = [1, 0.5, 0.15];
    
    %% Plot
    figH = figure(1);
    PlotYyWrapper(figH, TInfo, PrecipInfo, EmissionPotentialInfo, 'NorthEast');
    hgsave(sprintf('fig/oug_flux_c_rem_%d_days_lambda_%4.3f.fig', ...
        ceil(nT / TimeParams.intervalsPerDay), lambda));
%     print(figH, '-dpng', '-r0', sprintf('fig/oug_flux_c_rem_%d_days_lambda_%4.3f.png', ...
%         ceil(nT / TimeParams.intervalsPerDay), lambda));

    figH = figure(2);
    figPos = [100, 100, 500, 250];
    PlotDoubleWrapper(figH, TInfo, PrecipInfo, LeachateFluxInfo, 'NorthEast', figPos);
    % PrecipInfo - 'b'; LeachateFluxInfo - 'r'
    hgsave(sprintf('fig/precip_leachate_flux_%d_days_lambda_%4.3f.fig', ...
        ceil(nT / TimeParams.intervalsPerDay), lambda));
%     print(figH, '-dpng', '-r0', sprintf('fig/precip_leachate_flux_%d_days_lambda_%4.3f.png', ...
%         ceil(nT / TimeParams.intervalsPerDay), lambda));
    
    figH = figure(3);
    figPos = [200, 200, 500, 300];
    PlotYyWrapper(figH, TInfo, LeachateFluxInfo, LeachateConcentrationInfo, 'NorthEast', figPos);
    hgsave(sprintf('fig/concentr_leachate_flux_%d_days_lambda_%4.3f.fig', ...
        ceil(nT / TimeParams.intervalsPerDay), lambda));
%     print(figH, '-dpng', '-r0', sprintf('fig/concentr_leachate_flux_%d_days_lambda_%4.3f.png', ...
%         ceil(nT / TimeParams.intervalsPerDay), lambda));
    
    figH = figure(4);
    figPos = [300, 300, 500, 300];
    PlotDoubleWrapper(figH, TInfo, LeachateConcentrationInfo, EmissionPotentialInfo, ...
        'NorthEast', figPos);
    % LeachateConcentrationInfo - 'r', EmissionPotentialInfo - 'def'
    hgsave(sprintf('fig/concentr_c_rem_%d_days_lambda_%4.3f.fig', ...
        ceil(nT / TimeParams.intervalsPerDay), lambda));
%     print(figH, '-dpng', '-r0', sprintf('fig/concentr_c_rem_%d_days_lambda_%4.3f.png', ...
%         ceil(nT / TimeParams.intervalsPerDay), lambda));

    return
    
    %% Some generic plotting functions
    function PlotYyWrapper(figH, X, Var1, Var2, legendLocation, figPos)
        if nargin > 5
            set(figH, 'Position', figPos);
        end
        [axH, lH1, lH2] = plotyy(X.data, Var1.data, X.data, Var2.data);
        legend({Var1.name, Var2.name}, 'Location', legendLocation);
        xlabel(X.axisLabel);
        set(get(axH(1), 'ylabel'), 'string', Var1.axisLabel);
        set(get(axH(2), 'ylabel'), 'string', Var2.axisLabel);
        if isfield(Var1, 'color')
            set(lH1, 'color', Var1.color);
        end
        if isfield(Var2, 'color')
            set(lH2, 'color', Var2.color);
        end
    end

    function PlotDoubleWrapper(figH, X, Var1, Var2, legendLocation, figPos)
        if nargin > 5
            set(figH, 'Position', figPos);
        end
        if isfield(Var1, 'color')
            plot(t, Var1.data, 'Color', Var1.color);
        else
            plot(t, Var1.data);
        end
        hold on;
        if isfield(Var2, 'color')
            plot(t, Var2.data, 'Color', Var2.color);
        else
            plot(t, Var2.data);
        end
        hold off;
        xlabel(X.axisLabel);
        ylabel(Var1.axisLabel)
        legend({Var1.name, Var2.name}, 'Location', legendLocation);
    end
end