function CheckResults(ModelOutput, action, FILE_NAME, COMP_VARS)

    global EPSILON NUM_SIGMAS 
    global NO_VALIDATION SAVE_RESULTS COMPARE_RESULTS

    % Unpack some fields
    nT = ModelOutput.nT;
    
    % Remaining emission potential
    mRemRes = sum(ModelOutput.mRemaining(:, 2:nT+1), 1);
    ModelOutput.mRemRes = mRemRes;

    % Error check
    if ~RealEq(sum(ModelOutput.mIni) - sum(ModelOutput.mOutTotal), mRemRes(end), EPSILON)
        warning('ResultCheck:MassBalanceError', 'Absolute error is too high: err = %3.2e', ...
            abs(abs(sum(ModelOutput.mIni) - sum(ModelOutput.mOutTotal) - mRemRes(end))));
    end

    % Validate against previous runs
    if (action == SAVE_RESULTS)
        varList = {'t', 'cOutTotal', 'mOutTotal', 'qOutTotal', 'mRemRes'};
        for var = varList
            varName = var{:};
            eval(sprintf('%s = ModelOutput.%s;', varName, varName));
        end
        save(FILE_NAME, varList{:});
    elseif (action == COMPARE_RESULTS)
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