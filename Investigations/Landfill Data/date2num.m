function t = date2num(cellArr)
    [t_space{1:length(cellArr(:,2)),1}] = deal(' ');
    t = datenum([char(cellArr(:,2)) char(t_space) char(cellArr(:,3))]);
end