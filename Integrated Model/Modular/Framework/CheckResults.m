function CheckResults(ModelOutput, ModelParams, action, FILE_NAME, COMP_VARS)

    Const = DefineConstants();

    iSolute = 24;
    
    % Unpack some fields
    nT = ModelOutput.nT;
    
    % Remaining emission potential
    mRemRes = sum(ModelOutput.mRemaining(:, 2:nT+1, iSolute), 1);
    ModelOutput.mRemRes = mRemRes(:, :);

    mIniSum = sum(ModelOutput.mIni(:, :, iSolute), 1);
    mOutSum = cumsum(ModelOutput.mOutTotal(:, :, iSolute), 2);
    % Error check
    iErr = find(~RealEq(mIniSum - mOutSum, mRemRes, Const.EPSILON), 1, 'first');
    if ~isempty(iErr)
        warning('ResultCheck:MassBalanceError', ['Absolute mass balance error is too high ' ...
            'starting at step %d, solute index %d'], iErr, iSolute);
    else
        fprintf('Congratulations, mass balance is good!\n');
    end

    % Validate against previous runs
    if (action == Const.SAVE_RESULTS)
        varList = {'t', 'nT', 'mIni', 'qIn', 'cOutTotal', 'mOutTotal', 'qOutTotal', 'mRemRes', ...
            'mRemaining', 'cRemaining', 'cAll'};
        for var = varList
            varName = var{:};
            eval(sprintf('%s = ModelOutput.%s;', varName, varName));
        end
        varList = {varList{:}, 'ModelParams'};
        save(FILE_NAME, varList{:});
    elseif (action == Const.COMPARE_RESULTS)
        BaselineRes = load(FILE_NAME);
        nT = min(size(ModelOutput.cOutTotal, 2), size(BaselineRes.cOutTotal, 2));
        for var = COMP_VARS
            varName = var{:};
            DiffBl.(varName) = ModelOutput.(varName)(:, 1:nT, :) - ...
                BaselineRes.(varName)(:, 1:nT, :);
        end
        fprintf('Error analysis:\n');
        for varIdx = 1:numel(COMP_VARS)
            var = COMP_VARS{varIdx};
            fprintf('\t%s : %f\n', var, max(max(abs(DiffBl.(var)))));
        end
    end

end