function ftot = inhibition(Rp, Comp, CT, C)
nreac = length(Rp.max);

% Substrate inhibition factor for methanogenesis and sulfate reduction
Ks = Rp.inhib(1:nreac,1:2); fs = ones(nreac,1);
fs(2) = CT(2)/(CT(2) + Ks(2,1)); % total acid concentration 
fs(3) = CT(8)/(CT(8) + Ks(3,2)); % total acid concentration 

% pH inhibition factor for hydrolysis and methanogenesis
H_con = C(find(strcmp('H+.con',Comp.out))); fpH = ones(nreac,1);   Ki_pH = Rp.inhib(1:nreac,3);
k1 = find(Ki_pH ~= 0); 
fpH(k1) = Ki_pH(k1).^2./(Ki_pH(k1).^2+H_con^2);    % Siegrist et al 2002

% Total VFA inhibition for hydrolysis
fVFA = ones(nreac,1);   Ki_VFA = Rp.inhib(1:nreac,4); k2 = find(strcmp('H[Acetate].tot',Comp.master)); 
k1 = find(Ki_VFA ~= 0); 
fVFA(k1) = Ki_VFA(k1)./(CT(k2)+Ki_VFA(k1)); % Siegrist et al 2002

% Ammonia inhibition factor for methanogenesis
fT = ones(nreac,1); Ki_NH3 = Rp.inhib(1:nreac,5); k2 = find(strcmp('NH3.con',Comp.out));
k1 = find(Ki_NH3 ~= 0);
fT(k1) = Ki_NH3(k1).^2./(Ki_NH3(k1).^2+C(k2)^2);  % Siegrist et al 2002

% Sulfide inhibition factor for methanogenesis and sulfate reduction
fhs = ones(nreac,1); Ki_H2S = Rp.inhib(1:nreac,6); k2 = find(strcmp('H2S.con',Comp.out));
k1 = find(Ki_H2S ~= 0);
fhs(k1) = Ki_H2S(k1)./(Ki_H2S(k1)+C(k2));
k3 = find(strcmp('H[Acetate].con',Comp.out));
fhs(3) = (1-C(k2)/(547/(34.08*1000)))^0.401*((1+C(k3)/(54/(60.03*1000)))^1.08)^-1;

% % try out pH inhib for ammonia
% H_con = C(find(strcmp('H+.con',Comp.out))); Ki_pH = 1e-6;
% fpH(7) = Ki_pH.^2./(Ki_pH.^2+H_con^2);    % Siegrist et al 2002
% 
% % try out mass transfer inhibition for methanogenesis
% fm = ones(nreac,1); kmax = 1;
% if CT(2) > kmax
%     fm(2) = 1;
% else
%     fm(2) = (CT(2)/kmax);
% %     fm(2) = 0.01;
% end
% Total inhibition factor for each reaction
ftot = fs.*fpH.*fT.*fVFA.*fhs;
% ftot = fs.*fpH.*fT.*fVFA.*fhs.*fm;
end