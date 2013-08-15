function CheckResults(ModelOutput, action, FILE_NAME, COMP_VARS)

    Const = DefineConstants();

    % Unpack some fields
    nT = ModelOutput.nT;
    
    % Remaining emission potential
    mRemRes = sum(ModelOutput.mRemaining(:, 2:nT+1, 1), 1);
    ModelOutput.mRemRes = mRemRes(:, :, 1);

    mIniSum = sum(ModelOutput.mIni(:, :, 1));
    mOutSum = sum(ModelOutput.mOutTotal(:, :, 1));
    % Error check
    if ~RealEq(mIniSum - mOutSum, mRemRes(end), Const.EPSILON)
        warning('ResultCheck:MassBalanceError', 'Absolute error is too high: err = %3.2e', ...
            abs(abs(mIniSum - mOutSum - mRemRes(end))));
    end

    % Validate against previous runs
    if (action == Const.SAVE_RESULTS)
        varList = {'t', 'cOutTotal', 'mOutTotal', 'qOutTotal', 'mRemRes', 'cAll'};
        for var = varList
            varName = var{:};
            eval(sprintf('%s = ModelOutput.%s;', varName, varName));
        end
        save(FILE_NAME, varList{:});
    elseif (action == Const.COMPARE_RESULTS)
        BaselineRes = load(FILE_NAME);
        nEl = min(numel(ModelOutput.cOutTotal), numel(BaselineRes.cOutTotal));
        for var = COMP_VARS
            varName = var{:};
            DiffBl.(varName) = ModelOutput.(varName)(1:nEl) - BaselineRes.(varName)(1:nEl);
        end
        fprintf('Error analysis:\n');
        for varIdx = 1:numel(COMP_VARS)
            var = COMP_VARS{varIdx};
            fprintf('\t%s : %f\n', var, max(abs(DiffBl.(var))));
        end
    end

end