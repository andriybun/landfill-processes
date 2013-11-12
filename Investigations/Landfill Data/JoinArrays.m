function res = JoinArrays(arr1, colIdx1, arr2, colIdx2)
    % Function performs a left join of two tables arr1 and arr2. The tables are combined using 
    % fields colIdx1 and colIdx2 respectively.
    
    % Get sizes of tables
    [r1, c1] = size(arr1.data);
    [r2, c2] = size(arr2.data);
    
    % Column indices from table 2 to be copied into resulting table (all indices except colIdx2)
    iCopyCols2 = 1:c2;
    iCopyCols2(colIdx2) = [];
    
    % Initialize resulting table
    res = struct();
    res.data = nan(r1, c1+c2-1);
    res.data(:, 1:c1) = arr1.data;
    
    % Main loop
    for j = 1:r2
        % Find index of a row (i) of occurrence of j-th element from table 2 in table 1
        i = find(arr1.data(:, colIdx1) == arr2.data(j, colIdx2), true);
        if ~isempty(i)
            % Assign corresponding values from table 2 into resulting table
            res.data(i, c1+1:end) = arr2.data(j, iCopyCols2);
        end
    end

    % Truncate array - cut off NaN's at the beginning
    firstNumber = r1;
    for j = (c1+1):(c1+c2-1)
        firstNumber = min(firstNumber, find(~isnan(res.data(:, j)), true, 'first'));
    end
    res.data(1:firstNumber-1, :) = [];
    
    % Join headers
    res.headers = cat(1, arr1.headers, arr2.headers(iCopyCols2));
end