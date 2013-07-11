function CheckResults(ModelOutput, action, BASELINE_FILE_NAME, COMP_VARS)

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
        varList = {'cOutTotal', 'mOutTotal', 'qOutTotal', 'mRemRes'};
        for var = varList
            eval(sprintf('%s = ModelOutput.%s;', var{:}, var{:}));
        end
        save(BASELINE_FILE_NAME, varList{:});
    elseif (action == COMPARE_RESULTS)
        BaselineRes = load(BASELINE_FILE_NAME);
        nEl = min(numel(ModelOutput.cOutRes), numel(BaselineRes.cOutRes));
        DiffBl.cOutRes = ModelOutput.cOutRes(1:nEl) - BaselineRes.cOutRes(1:nEl);
        DiffBl.mOutRes = ModelOutput.mOutTotal(1:nEl) - BaselineRes.mOutRes(1:nEl);
        DiffBl.mRemRes = ModelOutput.mRemRes(1:nEl) - BaselineRes.mRemRes(1:nEl);
        
        fprintf('Error analysis:\n');
        for varIdx = 1:numel(COMP_VARS)
            var = COMP_VARS{varIdx};
            fprintf('\t%s : %f\n', var, max(abs(DiffBl.(var))));
        end
    end

end