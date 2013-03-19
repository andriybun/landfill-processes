function CheckPhaseTransport(ModelDim, SoilPar)
    % Perform consistency check of inputs
    errorCode = 0;

    if (size(SoilPar.d, 2) ~= SoilPar.nSolutes)
        % Incorrect definition of diffusion coefficients
        errorCode = 101;
    end

    if (ModelDim.nPhases ~= 2)
        % Wrong number of phases
        errorCode = 102;
    end

    if errorCode
        error('Error #%d\n', errorCode);
    end
end