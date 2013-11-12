function [out] = Store_Orchestra_Results(~,~,flag)

global MS tS timer_0 timer_flag
persistent MS_o tS_o 

if strcmp(flag,'init') == 1
    MS_o = []; tS_o = []; timer_0 = tic;
end
    
if strcmp(flag,'done')==0 && strcmp(flag,'init')==0
    MS_o = [MS_o;MS'];
    tS_o = [tS_o tS];
    out = timer_flag;
end

if strcmp(flag,'done') == 1
    save('Results_Orchestra','MS_o','tS_o')
end
end

