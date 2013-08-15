clear all

parameterName = 'lambda'; % 'beta', 'kExch', 'kExchPart', 'lambda'
outputParameterName = 'cOutTotal';
% outputParameterName = 'mRemRes';

FILENAME_TEMPLATE = '../Data/results_%s_%s.mat';

fileNameFilter = sprintf(FILENAME_TEMPLATE, parameterName, '*');
dirName = fileparts(FILENAME_TEMPLATE);

fileList = dir(fileNameFilter);
nScenarios = numel(fileList);
parameterValueListDouble = nan(1, nScenarios);
parameterValueListStr = cell(1, nScenarios);

for iScenario = 1:nScenarios
    fileName = fileList(iScenario).name;
    [~, fileName, ~] = fileparts(fileName);
    fileNameParts = strsplit(fileName, '_');
    parameterValueListStr{iScenario} = fileNameParts{end};
    parameterValueListDouble(iScenario) = str2double(parameterValueListStr(iScenario));
end

[~, sOrder] = sort(parameterValueListDouble);

legendInfo = cell(nScenarios, 1);
colorPalette = zeros(nScenarios, 3);
colorPalette(:, 1) = 1;
colorPalette(:, 2) = linspace(0.2, 0.95, nScenarios);

figure(1);

close all

hold on;

for iScenario = 1:nScenarios
    fileName = fileList(sOrder(iScenario)).name;
    
    data = load(fileName);
    
    nEl = numel(data.(outputParameterName));
    
    plot(data.t, data.(outputParameterName), 'Color', colorPalette(iScenario, :));

    legendInfo{iScenario} = sprintf('Sc. #%d: %s = %s', ...
        iScenario, parameterName, parameterValueListStr{sOrder(iScenario)});
end

ylabel('concentration [m^3/m^3]');
xlabel('time [days]');

hold off;

legend(legendInfo);