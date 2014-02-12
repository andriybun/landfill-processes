close all

% Starting date of dataset
sdate = '2012-01-01';
% sdate = '2012-05-16';

% SQL query for extracting data from database
QueryTmpl = [ ...
    'SELECT value, date, time\n' ...
    'FROM analysis\n' ...
    'WHERE methodno = %d AND measurement_pointno = %d AND date > "' sdate '"\n' ...
    'ORDER BY date, time'];

% Point index from 1 to 4. They are located in different landfills (5 - Braambergen, 6 -
% Wieringermeer). Mapping of points to landfill (to get weather data)
pointMap = [6, 5, 5, 5];
point1 = 3;
point2 = 4;

locationNameMap = {'Braambergen', 'Wieringermeer'};
locationName = locationNameMap{pointMap(point1) - 4};

% KNMI data
meteoDataFile = sprintf('data/KNMI_%s.csv', locationName);
knmiRaw = dlmread(meteoDataFile, ';');
knmiRaw(:, 2) = datenum(int2str(knmiRaw(:, 2)),'yyyymmdd');
% Get index of first day of interrest
is = find(knmiRaw(:, 2) == datenum(sdate), 1);
% Cut off all data before the start day and unnecessary columns:
%   2 - date
%   23 - precip (0.1 mm)
%   41 - evap (0.1 mm))
knmi = struct();
knmi.headers = {'time'; 'precipitation'; 'evaporation'; ...
    'net infiltration'; 'cumulative precip.'; 'cumulative evap.'; 'cumulative net inf.'};
knmi.data = knmiRaw(is:end, [2, 23, 41]);
netInf = NetInfiltration(knmi.data(:, 2), knmi.data(:, 3));
knmi.data = [knmi.data, netInf, cumsum([knmi.data(:, 2:3), netInf])];
knmi.data(:, 2:end) = knmi.data(:, 2:end) / 10;
save(sprintf('data/%sMeteoDailyProcessed', locationName), '-struct', 'knmi');

% Actual data
hoursRaw = fetch(dbConn, sprintf(QueryTmpl, 2, point1));
hoursRaw2 = fetch(dbConn, sprintf(QueryTmpl, 2, point2));
levelRaw = fetch(dbConn, sprintf(QueryTmpl, 3, point1));
levelRaw2 = fetch(dbConn, sprintf(QueryTmpl, 3, point2));
ecRaw = fetch(dbConn, sprintf(QueryTmpl, 5, point1));
ecRaw2 = fetch(dbConn, sprintf(QueryTmpl, 5, point2));
pumpRaw = fetch(dbConn, sprintf(QueryTmpl, 6, point1));
pumpRaw2 = fetch(dbConn, sprintf(QueryTmpl, 6, point2));
rainRaw = fetch(dbConn, sprintf(QueryTmpl, 7, pointMap(point1)));
rainRaw2 = fetch(dbConn, sprintf(QueryTmpl, 7, pointMap(point2)));

% Convert to numeric format
rain.data = [date2num(rainRaw), cell2mat(rainRaw(:, 1))];
rain.headers = {'time'; 'rain'};
level.data = [date2num(levelRaw), cell2mat(levelRaw(:, 1))];
level.headers = {'time'; 'level'};
level2.data = [date2num(levelRaw2), cell2mat(levelRaw2(:, 1))];
level2.headers = {'time'; 'level 2'};
hours.data = [date2num(hoursRaw), cell2mat(hoursRaw(:, 1))];
hours.headers = {'time'; 'hours'};
hours2.data = [date2num(hoursRaw2), cell2mat(hoursRaw2(:, 1))];
hours2.headers = {'time'; 'hours 2'};
pump.data = [date2num(pumpRaw), cell2mat(pumpRaw(:, 1))];
pump.headers = {'time'; 'pump'};
pump2.data = [date2num(pumpRaw2), cell2mat(pumpRaw2(:, 1))];
pump2.headers = {'time'; 'pump 2'};
ec.data = [date2num(ecRaw), cell2mat(ecRaw(:, 1))];
ec.headers = {'time'; 'ec'};
ec2.data = [date2num(ecRaw2), cell2mat(ecRaw2(:, 1))];
ec2.headers = {'time'; 'ec 2'};

% % Plotting
% figure(1)
% plotyy(rain.data(:, 1), rain.data(:, 2), pump.data(:, 1), pump.data(:, 2));
% figure(2)
% plotyy(level.data(:, 1), level.data(:, 2), pump.data(:, 1), pump.data(:, 2));

joined = JoinArrays(rain, 1, hours, 1);
joined = JoinArrays(joined, 1, hours2, 1);
joined = JoinArrays(joined, 1, pump, 1);
joined = JoinArrays(joined, 1, pump2, 1);
joined = JoinArrays(joined, 1, level, 1);
joined = JoinArrays(joined, 1, level2, 1);

LeachData = struct();
LeachData.headers = {'time', 'rainfall', 'pump hours 1', 'pump hours 2', ...
    'total debit 1', 'total debit 2'}';
LeachData.data = nan(size(joined.data, 1), 6);
LeachData.data(:, 1:6) = joined.data(:, 1:6);
save(sprintf('data/%sProcessed', locationName), '-struct', 'LeachData');

LevelData = struct();
LevelData.data = nan(size(joined.data, 1), 3);
LevelData.data(:, 1) = joined.data(:, 1);
LevelData.data(:, 2) = joined.data(:, 7);
LevelData.data(:, 3) = joined.data(:, 8);
LevelData.headers = {'time', 'level 1', 'level 2'}';
save(sprintf('data/%sLevelProcessed', locationName), '-struct', 'LevelData');

% save('data/RaindataProcessed', '-struct', 'RainData');

% % Joins
% JointQueryTmpl = [ ...
%     'SELECT *\n' ...
%     'FROM (%s) t1\nFULL JOIN (%s) t2\n' ...
%     'ON ((t1.date = t2.date) AND (t1.time = t2.time))\n' ...
%     'ORDER BY date, time'];
% JointQuery = sprintf(JointQueryTmpl, sprintf(QueryTmpl, 3, point), sprintf(QueryTmpl, 6, point));
% data = fetch(dbConn, JointQuery);
% data(1:10, :)
