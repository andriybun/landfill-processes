function LearnGui
    addpath('../../Common');

    params = struct();
    params.name = 'Just a name';
    params.date = '25-07-2014';
    params.int = int32(365);
    params.pi = 3.1415;
    params.exp = exp(1);
    params.one = int32(1);
    params.x = 0:0.1:30;
    params.minVal = 10;
    params.maxVal = 20;
    
    strct = struct();
    strct.a = 1;
    strct.b = 2;
    params.struct = strct;
    
    runParams = struct();
    runParams.dir = cd;
    runParams.file = mfilename;
        
    result = RunWithGui('Test', @Dummy, params, strct, runParams);
    
    return
    
    function result = Dummy(input, strct, ~, prBar)
        result = struct();
        result.x = input.x;
        result.sin = sin(result.x);
        result.cos = cos(result.x);
        result.tan = tan(result.x);
        prBar.pvalue = 100;
        disp(input.struct);
    end
end