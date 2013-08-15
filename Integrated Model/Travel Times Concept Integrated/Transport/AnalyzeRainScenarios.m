clear all

parameterName = 'Rain_Data'; % 'beta', 'kExch', 'kExchPart', 'lambda'
% outputParameterName = 'cOutTotal';
outputParameterName = 'mRemRes';

FILENAME_TEMPLATE = '../Data/results_%s_%s.mat';

fileNameFilter = sprintf(FILENAME_TEMPLATE, parameterName, '*');
dirName = fileparts(FILENAME_TEMPLATE);

fileList = dir(fileNameFilter);
nScenarios = numel(fileList);
parameterValueListDouble = nan(1, nScenarios);
parameterNameListStr = cell(1, nScenarios);

for iScenario = 1:nScenarios
    fileName = fileList(iScenario).name;
    [~, fileName, ~] = fileparts(fileName);
    fileNameParts = strsplit(fileName, '_');
    scenarioNameWordList = fileNameParts(4:end);
    scenarioName = scenarioNameWordList{1};
    for iWord = 2:numel(scenarioNameWordList)
        word = scenarioNameWordList{iWord};
        scenarioName = [scenarioName ' ' word];
    end
    parameterNameListStr{iScenario} = scenarioName;
end

legendInfo = cell(nScenarios, 1);
colorPalette = zeros(nScenarios, 3);
colorPalette(:, 1) = 1;
colorPalette(:, 2) = linspace(0.2, 0.95, nScenarios);

figure(1);

close all

hold on;

for iScenario = 1:nScenarios
    fileName = fileList(iScenario).name;
    
    data = load(fileName);
    
    nEl = numel(data.(outputParameterName));
    
    plot(data.t, data.(outputParameterName), 'Color', colorPalette(iScenario, :));

    legendInfo{iScenario} = sprintf('Sc. #%d: %s = %s', ...
        iScenario, strrep(parameterName, '_', ' '), parameterNameListStr{iScenario});
end

ylabel('concentration [m^3/m^3]');
xlabel('time [days]');

hold off;

legend(legendInfo);