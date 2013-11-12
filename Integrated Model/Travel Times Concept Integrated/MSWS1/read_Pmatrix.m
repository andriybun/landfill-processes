function [MT_l, MT_g, OUTPUT, PARAMETERS, INHIBITION_TYPE] = read_Pmatrix(fn)
 
ii = 0;
fid = fopen(fn); tline = fgetl(fid); %Pmatrix = [];
while tline ~= -1,
    ii = ii+1;
    line = textscan(tline,'%s','delimiter',';');
    Pmatrix{ii,1} = line{:}';
    tline = fgetl(fid);
end
fclose (fid);

for j = 1:length(Pmatrix(:))
    
    if isempty(Pmatrix{j}{1}) == 1
        break
    end
    
    for i = 1:length(Pmatrix{j})
        if isempty(Pmatrix{j}{i}) == 1
            break
        else
            MT_l{j,i} = Pmatrix{j}{i};
        end
    end
    N = j;
end

ii = 0; 
for j = N+1:length(Pmatrix(:))
    
    if isempty(Pmatrix{j}{1}) == 1 && ii > 0
        break
    elseif isempty(Pmatrix{j}{1}) == 1 
        N = N+1;
    end
    
    if isempty(Pmatrix{j}{1}) == 0
        ii = 1;
        for i = 1:length(Pmatrix{j})
            if isempty(Pmatrix{j}{i}) == 1
                break
            else
                OUTPUT{j-N,i} = Pmatrix{j}{i};
            end
        end
    end
end

ii = 0; W = length(OUTPUT(1,:));
for j = N+1:length(Pmatrix(:))
    
    for i = 1+W:length(Pmatrix{j})
        
        if isempty(Pmatrix{j}{i}) == 0
            ii = 1;
            PARAMETERS{j-N,i-W} = Pmatrix{j}{i};
        end
        
        if isempty(Pmatrix{j}{i}) == 1 && ii == 0
            W = W+1;
        end
        
        if isempty(Pmatrix{j}{i}) == 1 && ii == 1
            break
        end
    end
end

ii = 0; W = W + length(PARAMETERS(1,:));
for j = N+1:length(Pmatrix(:))
    
    for i = 1+W:length(Pmatrix{j})
        
        if isempty(Pmatrix{j}{i}) == 0
            ii = 1;
            INHIBITION_TYPE{j-N,i-W} = Pmatrix{j}{i};
        end
        
        if isempty(Pmatrix{j}{i}) == 1 && ii == 0
            W = W+1;
        end
        
        if isempty(Pmatrix{j}{i}) == 1 && ii == 1
            break
        end
    end
end

ii = 0; W = W + length(INHIBITION_TYPE(1,:));
for j = N+1:length(Pmatrix(:))
    
    for i = 1+W:length(Pmatrix{j})
        
        if isempty(Pmatrix{j}{i}) == 0
            ii = 1;
            MT_g{j-N,i-W} = Pmatrix{j}{i};
        end
        
        if isempty(Pmatrix{j}{i}) == 1 && ii == 0
            W = W+1;
        end
        
        if isempty(Pmatrix{j}{i}) == 1 && ii == 1
            break
        end
    end
end   
end
