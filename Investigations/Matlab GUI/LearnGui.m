function LearnGui
    addpath('../../Common');

    input = struct();
    input.name = 'Just a name';
    input.date = '25-07-2014';
    input.int = int32(365);
    input.pi = 3.1415;
    input.exp = exp(1);
    input.one = int32(1);
    input.x = 0:0.1:30;
    input.minVal = 10;
    input.maxVal = 20;
    
    strct = struct();
    strct.a = 1;
    strct.b = 2;
    input.struct = strct;
    
    runParams = struct();
    runParams.dir = cd;
    runParams.file = mfilename;
        
    % Here 3 parameter and on are passed to a function specified by a handle (@Dummy in this case)
    RunWithGui('Test', @Dummy, input, strct, runParams);
    
    return
    
    function result = Dummy(input, strct, runParams, prBar)
        % This function is executed by pressing button "Start". All parameters except last one are
        % the same as passed above. The last one is progress bar. Can be set by changing value of
        % pvalue property within range 0..100
        result = struct();
        result.x = input.x;
        result.sin = sin(result.x);
        result.cos = cos(result.x);
        result.tan = tan(result.x);
        % To illustrate plots of slices of 3d arrays or 3d meshes
        result.all = permute([result.sin; result.cos; result.tan], [3, 2, 1]);
        prBar.pvalue = 100;
        disp(input.struct);
    end
end