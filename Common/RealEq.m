function res = RealEq(val1, val2, EPSILON)
    res = abs(val1 - val2) < EPSILON;
end