function ShowPlots(qOutTotal, mOutTotal, emissionPotential, rainData, lambda, TimeParams)

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
    EmissionPotentialInfo.data = emissionPotential;
    EmissionPotentialInfo.name = 'Emission potential';
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
    close all;
    PlotYyWrapper(TInfo, PrecipInfo, EmissionPotentialInfo, 'NorthEast');

    figPos = [100, 100, 500, 250];
    PlotDoubleWrapper(TInfo, PrecipInfo, LeachateFluxInfo, 'NorthEast', figPos);
    
    figPos = [200, 200, 500, 300];
    PlotYyWrapper(TInfo, LeachateFluxInfo, LeachateConcentrationInfo, 'NorthEast', figPos);
    
    figPos = [300, 300, 500, 300];
    PlotDoubleWrapper(TInfo, LeachateConcentrationInfo, EmissionPotentialInfo, ...
        'NorthEast', figPos);

    return
    
    %% Some generic plotting functions
    function PlotYyWrapper(X, Var1, Var2, legendLocation, figPos)
        figH = figure();
        if nargin > 4
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
        % Save to file
        hgsave(figH, GenerateFileName(Var1, Var2, 'fig'));
%         print(figH, '-dpng', '-r0', GenerateFileName(Var1, Var2, 'png'));
    end

    function PlotDoubleWrapper(X, Var1, Var2, legendLocation, figPos)
        figH = figure();
        if nargin > 4
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
        % Save to file
        fileName = GenerateFileName(Var1, Var2, 'fig');
        hgsave(figH, fileName);
%         print(figH, '-dpng', '-r0', GenerateFileName(Var1, Var2, 'png'));
    end

    function fileName = GenerateFileName(Var1, Var2, ext)
        fileName = sprintf('fig/%s_-_%s_%d_days_lambda_%4.3f.%s', ...
            strrep(Var1.name, ' ', '_'), strrep(Var2.name, ' ', '_'), ...
            TimeParams.maxDays, lambda, ext);
    end
end