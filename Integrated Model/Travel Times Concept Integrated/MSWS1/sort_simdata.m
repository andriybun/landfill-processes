function data_sorted = sort_simdata(datas,ts,tm)

data_sorted = []; 
for i = 1:length(datas);
    
    data = datas{i}; t = ts{i}; tf = tm{i}; data_av = []; 
    for j = 1:length(tf)
        
        k1 = find(tf(j) == t);
        
        if k1 > 0
            data_av(j) = data(k1);
        else
            k1 = find(tf(j) > t);
        
            if k1 > 0
                k1 = k1(end);
                data_av(j) = (data(k1)+data(k1+1))/2;
            else
            data_av(j) = 0;
            end
        end
    end
    data_sorted = [data_sorted data_av];
end

