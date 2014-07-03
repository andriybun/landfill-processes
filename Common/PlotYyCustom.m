function figH = PlotYyCustom(X, Var1, Var2, legendLocation, figPos)
        figH = figure();
        if nargin > 4
            set(figH, 'Position', figPos);
        end
        nLines = size(Var2.data, 2);
        [axH, lH1, lH2] = plotyy(X.data, Var1.data, X.data, Var2.data);
        legend({Var1.name{1}, Var2.name{:}}, 'Location', legendLocation);
        xlabel(X.axisLabel);
        % set(axH, {'ycolor'}, {Var1.color; Var2.color(1, :)});
        set(axH, {'ycolor'}, {'black'; 'black'});
        set(get(axH(1), 'ylabel'), 'string', Var1.axisLabel);
        set(get(axH(2), 'ylabel'), 'string', Var2.axisLabel);
        if isfield(Var1, 'color')
            set(lH1, 'color', Var1.color);
            set(lH1, 'LineStyle', Var1.lineStyle);
        end
        if isfield(Var2, 'color')
            for iLine = 1:nLines
                set(lH2(iLine), 'color', Var2.color(iLine, :));
                set(lH2(iLine), 'LineStyle', Var2.lineStyle{iLine});
            end
        end
%         xlim(axH(1), [0, 200]);
%         xlim(axH(2), [0, 200]);
        ylim(axH(1), [0, 1.1 * max(Var1.data)]);
        ylim(axH(2), [0, 1.1 * max(max(Var2.data))]);
        % Save to file
%         hgsave(figH, GenerateFileName(Var1, Var2, 'fig'));
%         print(figH, '-dpng', '-r0', GenerateFileName(Var1, Var2, 'png'));
end