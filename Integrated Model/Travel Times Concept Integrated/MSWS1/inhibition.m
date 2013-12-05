function I = inhibition(Rp, Comp, MT, C, V)
nreac = length(Rp.max(:,1)); CT = MT/V; CT = [CT; C];    
f1 = ones(nreac,1); f2 = ones(nreac,1); f3 = ones(nreac,1); f4 = ones(nreac,1); f5 = ones(nreac,1);

% f1: Substrate inhibition 
k1 = find(Rp.inhib(:,1) == 1);     
Ks = Rp.inhib(k1,4); 
Ci = CT(Rp.inhib(k1,3));
for i = 1:length(k1)
    f1(Rp.inhib(k1(i),2)) = f1(Rp.inhib(k1(i),2))*(Ci(i)./(Ci(i) + Ks(i))); % total acid concentration
end

% f2:
k1 = find(Rp.inhib(:,1) == 2);     
Ki = Rp.inhib(k1,4); 
Ci = CT(Rp.inhib(k1,3));
f2(Rp.inhib(k1,2)) = Ki./(Ci + Ki); % Siegrist et al 2002
    
% f3:
k1 = find(Rp.inhib(:,1) == 3);     
Ki = Rp.inhib(k1,4); 
Ci = CT(Rp.inhib(k1,3));
f3(Rp.inhib(k1,2)) = Ki.^2./(Ki.^2 + Ci.^2);    % Siegrist et al 2002
    
% f4:
k1 = find(Rp.inhib(:,1) == 4);     
Ki = Rp.inhib(k1,4); 
Ci = CT(Rp.inhib(k1,3));
f4(Rp.inhib(k1,2)) = (1-Ci./Ki).^0.401;  

% f5:
k1 = find(Rp.inhib(:,1) == 5);     
Ki = Rp.inhib(k1,4); 
Ci = CT(Rp.inhib(k1,3));
f5(Rp.inhib(k1,2)) = (1+(Ci./Ki).^1.08).^-1;  

% I
I = f1.*f2.*f3.*f4.*f5;
end