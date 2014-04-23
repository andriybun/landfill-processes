function ShowPlots(ModelOutput, ModelParams, TimeParams, iSpecies, tShow)
 
    Const = DefineConstants();

    % Unpack resulting arrays
    nT = ModelOutput.nT;
    lambda = ModelParams.lambda;
    qOutTotal = ModelOutput.qOutTotal;
    mOutTotal = ModelOutput.mOutTotal;
    cOutTotal = ModelOutput.cOutTotal;
%     emissionPotential = sum(ModelOutput.mRemaining(:, 2:nT+1, iSpecies), 1);

    if nargin < 6
        nT = TimeParams.numIntervals;
        tShow = true(1, nT);
        t = TimeParams.t(1:nT);
    else
        t = TimeParams.t(tShow);
    end

    %% Prepare data structures to be plotted:
    TInfo = struct();
    TInfo.data = t;
    TInfo.axisLabel = 'time [days]';
    
    PrecipInfo = struct();
    PrecipInfo.data = ModelOutput.qIn(:, tShow);
    PrecipInfo.name = {'Precipitation'};
    PrecipInfo.axisLabel = 'precipitation [m/hour]';
    PrecipInfo.color = [0.4, 0.7, 1];
    
%     EmissionPotentialInfo = struct();
%     EmissionPotentialInfo.data = emissionPotential(:, tShow);
%     EmissionPotentialInfo.name = {'Emission potential'};
%     EmissionPotentialInfo.axisLabel = 'emission potential m^3]';
%     EmissionPotentialInfo.color = [0.4, 0.4, 0.4];
    
    LeachateFluxInfo = struct();
    LeachateFluxInfo.data = qOutTotal(1, tShow);
    LeachateFluxInfo.name = {'Leachate volume flux'};
    LeachateFluxInfo.axisLabel = 'out flux [m/hour]';
    LeachateFluxInfo.color = [0.55, 0.25, 0.08];
    
    ImmobileConcentrationInfo = struct();
    pv = [ModelParams.totalPv - ModelParams.totalPv / (1 + ModelParams.beta); ...
        ModelParams.totalPv / (1 + ModelParams.beta)];
    ImmobileConcentrationInfo.data = (pv(1) * ModelOutput.cRemaining(1, 2:nT+1, iSpecies) + ...
        pv(2) * ModelOutput.cRemaining(2, 2:nT+1, iSpecies)) / ModelParams.totalPv;
    ImmobileConcentrationInfo.name = {'Remaining concentration'};
    ImmobileConcentrationInfo.axisLabel = 'concentration [kg/m^3]';
    ImmobileConcentrationInfo.color = [0.4, 0.4, 0.4];
    
    LeachateConcentrationInfo = struct();
    LeachateConcentrationInfo.data = cOutTotal(:, tShow, iSpecies);
    LeachateConcentrationInfo.data = cat(1, LeachateConcentrationInfo.data, ...
        ImmobileConcentrationInfo.data);
    LeachateConcentrationInfo.name = {'Leachate concentration', 'Emission potential'};
    LeachateConcentrationInfo.axisLabel = 'concentration [kg/m^3]';
    LeachateConcentrationInfo.color = [1, 0.5, 0.15];
    LeachateConcentrationInfo.color = cat(1, LeachateConcentrationInfo.color, ...
        ImmobileConcentrationInfo.color);
    
    %% Plot
    close all;
    PlotYyWrapper(TInfo, PrecipInfo, ImmobileConcentrationInfo, 'NorthEast');

    figPos = [100, 100, 500, 250];
    PlotDoubleWrapper(TInfo, PrecipInfo, LeachateFluxInfo, 'NorthEast', figPos);
    
    figPos = [200, 200, 500, 300];
    PlotYyWrapper(TInfo, LeachateFluxInfo, LeachateConcentrationInfo, 'NorthEast', figPos);
    
%     figPos = [300, 300, 500, 300];
%     PlotYyWrapper(TInfo, LeachateConcentrationInfo, ImmobileConcentrationInfo, ...
%         'NorthEast', figPos);

    figPos = [400, 400, 500, 250];
    PlotDoubleWrapper(TInfo, ImmobileConcentrationInfo, LeachateConcentrationInfo, ...
        'NorthEast', figPos);
    
    return
    
    %% Some generic plotting functions
    function PlotYyWrapper(X, Var1, Var2, legendLocation, figPos)
        figH = figure();
        if nargin > 4
            set(figH, 'Position', figPos);
        end
        nLines = size(Var2.data, 1);
        [axH, lH1, lH2] = plotyy(X.data, Var1.data, X.data, Var2.data);
        legend({Var1.name{1}, Var2.name{:}}, 'Location', legendLocation);
        xlabel(X.axisLabel);
        set(axH, {'ycolor'}, {Var1.color; Var2.color(1, :)});
        set(get(axH(1), 'ylabel'), 'string', Var1.axisLabel);
        set(get(axH(2), 'ylabel'), 'string', Var2.axisLabel);
        if isfield(Var1, 'color')
            set(lH1, 'color', Var1.color);
        end
        if isfield(Var2, 'color')
            set(lH2(1), 'color', Var2.color(1, :));
            set(lH2(1), 'LineStyle', '--');
            if (nLines > 1)
                set(lH2(2), 'color', Var2.color(2, :));
                set(lH2(2), 'LineStyle', '-.');
            end
        end
%         xlim(axH(1), [0, 200]);
%         xlim(axH(2), [0, 200]);
        ylim(axH(1), [0, 1.1 * max(Var1.data)]);
        ylim(axH(2), [0, 1.1 * max(max(Var2.data))]);
        % Save to file
        hgsave(figH, GenerateFileName(Var1, Var2, 'fig'));
        print(figH, '-dpng', '-r0', GenerateFileName(Var1, Var2, 'png'));
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
            plot(t, Var2.data(1, :), 'Color', Var2.color(1, :));
        else
            plot(t, Var2.data(1, :));
        end
        hold off;
        xlabel(X.axisLabel);
        ylabel(Var1.axisLabel)
        legend({Var1.name{1}, Var2.name{1}}, 'Location', legendLocation);
%         xlim([0, 200]);
        % Save to file
        fileName = GenerateFileName(Var1, Var2, 'fig');
        hgsave(figH, fileName);
        print(figH, '-dpng', '-r0', GenerateFileName(Var1, Var2, 'png'));
    end

    function fileName = GenerateFileName(Var1, Var2, ext)
        fileName = sprintf('fig/%s_-_%s_sp_#%02d_%d_days_lambda_%4.3f.%s', ...
            strrep(Var1.name{1}, ' ', '_'), strrep(Var2.name{1}, ' ', '_'), ...
            iSpecies, TimeParams.maxDays, lambda, ext);
    end
end