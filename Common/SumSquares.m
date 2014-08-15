function err_ = SumSquares(v1, v2, k)
    if nargin < 3
        v2Cum = cumsum(v2);
        pos = find(v2Cum > 0.9, 1, 'first');
        [nRows, nCols] = size(v2);
        k = ones(nRows, nCols);
        k(pos:nRows) = 1 - v2Cum(1:(nRows-pos+1));
    end
    err_ = sum((k .* v1 - v2) .^ 2);
end