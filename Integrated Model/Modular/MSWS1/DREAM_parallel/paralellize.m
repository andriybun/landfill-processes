function varargout = paralellize(fHandle, loopArr, numWorkers, varargin)
    % Analyze arguments
    passVarargin = ~isempty(varargin);
    loopArrIsVector = isvector(loopArr);
    loopArrIsColumn = false;
    
    % Check number of output arguments of a function passed
    numArgOut = nargout(fHandle);
    
    
    if loopArrIsVector
        num = numel(loopArr);
        loopArrIsColumn = iscolumn(loopArr);
        if ~loopArrIsColumn
            loopArr = loopArr';
        end
    else
        num = size(loopArr, 1);
    end

%     % Checking if there are any workers still running
%     matlabpoolIsOpen = (matlabpool('size') > 0);
%     % If yes, close them
%     if matlabpoolIsOpen
%         matlabpool close
%     end
% 
%     % Start matlabpool
%     matlabpool(min(numWorkers, num));

    % Compute workload per worker and start of batch for each worker
    numElemPerWorker = ceil(num / numWorkers);
    firstElement = 1:numElemPerWorker:num;

    % Parallel loop
    spmd
        % First and last elements of batch for current worker
        start = firstElement(labindex);
        finish = min(num, firstElement(labindex) + numElemPerWorker - 1);
        % Initialize cell array for all outputs
        allOut = cell(finish - start + 1, numArgOut);
        count = 1;
        % Process all values from batch
        for idx = start:finish
            % fprintf('worker %d: idx = %d\n', labindex, idx);
            if passVarargin
                [allOut{count, :}] = fHandle(loopArr(idx, :), varargin{:});
            else
                [allOut{count, :}] = fHandle(loopArr(idx, :));
            end
            % outputRaw(count) = allOut{1};
            count = count + 1;
        end
    end

    % Assemble results
    varargout = cell(1, numArgOut);
    
    coreIdx = 1;
    args = allOut{coreIdx};
    for argIdx = 1:numArgOut
        varargout{argIdx} = [args{:, argIdx}];
    end
    for coreIdx = 2:numel(allOut)
        args = allOut{coreIdx};
        for argIdx = 1:numArgOut
            varargout{argIdx} = cat(2, varargout{argIdx}, [args{:, argIdx}]);
        end
    end
    
    % Transpose results to column vectors, if needed
    if loopArrIsVector && loopArrIsColumn
        for argIdx = 1:numArgOut
            varargout{argIdx} = varargout{argIdx}';
        end
    end
    
%     % Stop workers
%     matlabpool close
end