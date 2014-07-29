function LearnGui
    addpath('../../Common');

    params = struct();
    params.name = 'Just a name';
    params.date = '25-07-2014';
    params.int = int32(365);
    params.pi = 3.1415;
    params.exp = exp(1);
    params.one = int32(1);
    
    result = RunWithGui('Test', params, @Dummy);
    
    return
    
    function result = Dummy(input, prBar)
        result = struct();
        result.x = 0:0.1:10;
        result.sin = sin(result.x);
        result.cos = cos(result.x);
        result.tan = tan(result.x);
        prBar.pvalue = 100;
    end
end