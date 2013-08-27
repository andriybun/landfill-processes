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
        warning('ResultCheck:MassBalanceError', ['Absolute mass balance error is too high: ' ...
            'err = %3.2e'], abs(abs(mIniSum - mOutSum - mRemRes(end))));
    else
        fprintf('Congratulations, mass balance is good!\n');
    end

    % Validate against previous runs
    if (action == Const.SAVE_RESULTS)
        varList = {'t', 'nT', 'mIni', 'cOutTotal', 'mOutTotal', 'qOutTotal', 'mRemRes', ...
            'mRemaining', 'cRemaining', 'cAll'};
        for var = varList
            varName = var{:};
            eval(sprintf('%s = ModelOutput.%s;', varName, varName));
        end
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