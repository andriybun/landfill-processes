function ModelOutput = SummarizeResults(ModelOutputAll, ModelParams)
    nCells = numel(ModelOutputAll);
    ModelOutput = ModelOutputAll{1};
    
    iCell = 1;
    ModelOutput.cRemaining = ModelOutputAll{iCell}.cRemaining * ModelParams{iCell}.totalPv;
    ModelOutput.cOutTotal = ModelOutputAll{iCell}.cOutTotal * ModelParams{iCell}.totalPv;
    sumPv = ModelParams{iCell}.totalPv;
    
    for iCell = 2:nCells
        ModelOutput.mIni = ModelOutput.mIni + ModelOutputAll{iCell}.mIni;
        ModelOutput.qIn = ModelOutput.qIn + ModelOutputAll{iCell}.qIn;
        ModelOutput.qOutTotal = ModelOutput.qOutTotal + ModelOutputAll{iCell}.qOutTotal;
        ModelOutput.mOutTotal = ModelOutput.mOutTotal + ModelOutputAll{iCell}.mOutTotal;
        ModelOutput.cRemaining = ModelOutput.cRemaining + ...
            ModelOutputAll{iCell}.cRemaining * ModelParams{iCell}.totalPv;
        ModelOutput.mRemaining = ModelOutput.mRemaining + ModelOutputAll{iCell}.mRemaining;
        ModelOutput.cOutTotal = ModelOutput.cOutTotal + ...
            ModelOutputAll{iCell}.cOutTotal * ModelParams{iCell}.totalPv;
        sumPv = sumPv + ModelParams{iCell}.totalPv;
    end
    
    ModelOutput.cRemaining = ModelOutput.cRemaining / sumPv;
    ModelOutput.cOutTotal = ModelOutput.cRemaining / sumPv;
end