function PreProcessLandfillData

    %% Leachate data
    LeachData = load('data/WieringermeerPPNPPZ');
    LeachData.data = flipud(LeachData.data);

    isNumber = ~isnan(LeachData.data);
    
    LeachDataNew = struct();
    LeachDataNew.headers = {'time', 'rainfall', 'pump hours 1', 'pump hours 2', ...
        'total debit 1', 'total debit 2'}';
    LeachDataNew.data = nan(size(LeachData.data, 1), 6);
    LeachDataNew.data(:, 1) = LeachData.data(:, 1);
    LeachDataNew.data(:, 2) = LeachData.data(:, 6);
    LeachDataNew.data(isNumber(:, 2), 3) = LeachData.data(isNumber(:, 2), 2);
    LeachDataNew.data(isNumber(:, 5), 3) = LeachData.data(isNumber(:, 5), 5);
    LeachDataNew.data(isNumber(:, 3), 4) = LeachData.data(isNumber(:, 3), 3);
    LeachDataNew.data(isNumber(:, 4), 4) = LeachData.data(isNumber(:, 4), 4);
    LeachDataNew.data(isNumber(:, 7), 5) = LeachData.data(isNumber(:, 7), 7);
    LeachDataNew.data(isNumber(:, 10), 5) = LeachData.data(isNumber(:, 10), 10);
    LeachDataNew.data(isNumber(:, 8), 6) = LeachData.data(isNumber(:, 8), 8);
    LeachDataNew.data(isNumber(:, 9), 6) = LeachData.data(isNumber(:, 9), 9);
    LeachDataNew.data(:, 5) = LeachDataNew.data(:, 5) - LeachDataNew.data(1, 5);
    LeachDataNew.data(:, 6) = LeachDataNew.data(:, 6) - LeachDataNew.data(1, 6);
    
    save('data/WieringermeerProcessed', '-struct', 'LeachDataNew');
    
    %% Level data
    LevelData = load('data/WieringermeerLevel');
    LevelData.data = flipud(LevelData.dataL);
    LevelData.headers = LevelData.headersL;
    LevelData.data(:, 4:end) = [];
    LevelData.headers(4:end) = [];
    LevelData = rmfield(LevelData, {'dataL', 'headersL'});

    save('data/WieringermeerLevelProcessed', '-struct', 'LevelData');
    
    %% Rainfall data
    RainData = load('data/raindata');
    
    RainData.data = [RainData.rf_time, RainData.rf_mean, RainData.ev_mean];
    RainData.data(:, 4) = NetInfiltration(RainData.data(:, 2), RainData.data(:, 3));
    
    RainData = rmfield(RainData, {'rf_mean', 'rf_time', 'ev_mean'});
    RainData.headers = {'time'; 'rainfall'; 'evaporation'; 'net infiltration'};

    save('data/RaindataProcessed', '-struct', 'RainData');
    
    %%
end