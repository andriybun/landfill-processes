function PlotAdd(fig, dataX, dataY, legendAdd)
    figure(fig);
    hold on;
    plot(dataX, dataY, 'r');
    lH = findobj(fig, 'Type', 'axes', 'Tag', 'legend');
    legendEntries = get(lH, 'String');
    legend([legendEntries, legendAdd]);
    hold off;
end