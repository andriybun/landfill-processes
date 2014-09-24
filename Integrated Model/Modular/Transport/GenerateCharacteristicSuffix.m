function fileNameSuffix = GenerateCharacteristicSuffix(ModelParams, ParameterOfInterest)
    fileNameSuffix = sprintf('beta_%3.2f_kEx_%4.3f_(%s)', ...
        ModelParams.beta, ModelParams.kExch, ParameterOfInterest.name);
end