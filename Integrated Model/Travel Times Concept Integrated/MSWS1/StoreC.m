function [status] = StoreC(t,y,x)

global C Vv Call V2 tt R Rall 
if strcmp(x,'done')==0 && strcmp(x,'init')==0
     
    Call  = [Call; C'];
    V2 = [V2 Vv(end)];
    tt = [tt;t(end)];
    Rall = [Rall; R'];
end
status = 0;
end

