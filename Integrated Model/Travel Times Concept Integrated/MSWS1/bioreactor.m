function dMTdt = bioreactor(t, MT, Comp, Pm, S, Rp, ORI, H)

[nreac,ncomp] = size(S); ng = 6; GT = MT(ncomp+1:ncomp+ng); V = MT(ncomp+ng+1); MT = MT(1:ncomp);
%---------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------
% MEASURES TO STABILIZE CALCULATIONS

% Kick out negative concentrations for ORCHESTRA calculations
MT_o = MT;
for i=1:ncomp
    if MT_o(i) < 1e-10 
        if i ~= find(strcmp('H+.tot',Comp.master))
           MT_o(i) = 1e-10;
        end
    end
end

% Make sure Biomass concentrations go to zero
for i = 10:12
    if MT(i) <= 1e-10
        MT(i) = 0;
    end
end
%---------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------
% CALCULATE SPECIATIONS OF MASTER SPECIES --> ORCHESTRA

global C Vv R
% n = (1:(ncomp+length(Comp.consti)));
% CT = [MT_o/V;Comp.consti'];
n = (1:ncomp);
CT = [MT_o/V];
ro = ORI.Calculate(n, CT);
C = ro(ncomp+length(Comp.consti)+length(Comp.inii)+1:end);
Vv = V;
%---------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------
% CALCULATE TOTAL INHIBITION FACTOR

ftot = inhibition(Rp, Comp, CT, C);
%---------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------
% DEFINE DERIVATIVES IN LIQUID PHASE

% Rate constant at t per compound per biochemical reaction (mol/mol/d)
%                           --> B0 = Stoichiometry x Rate constant x inhibition factor
B0 = (S.*(Rp.max.*ftot*ones(1,ncomp)));   

% Rate constant at t per compound per biochemical reaction times Ms or Mx (mol/d)
%                           --> B1 = B0 x Ms/Mx 
for i = 1:nreac;
    pC_dom = [1 11 12 10 11 12 4];
    for j = 1:ncomp;
             B1(i,j) = B0(i,j)*MT(pC_dom(i));
    end
end
dMTdt = sum(B1)'; % Total derivative per compound in liquid phase (all biochemical reactions summed)

%---------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------
%  CORRECT DERIVATIVES IN LIQUID PHASE FOR MASS TRANSFER

T = Pm(1); p = Pm(2); R = Pm(3); Vg = Pm(5); kla = Pm(6); rN2in = Pm(7)*1e-3*p/R/T; 

k1 = find(strcmp('H2CO3.con',Comp.out));    k2 = find(strcmp('H2CO3.tot',Comp.master));
k3 = find(strcmp('CH4',Comp.master));       k4 = find(strcmp('H2S.con',Comp.out));
k5 = find(strcmp('H2S.tot',Comp.master));   k6 = find(strcmp('NH3.con',Comp.out));
k7 = find(strcmp('NH3.tot',Comp.master));   k8 = find(strcmp('H2O.con',Comp.out));
k9 = find(strcmp('H2O',Comp.master));       

kk = [0 k1 k3 k4 k6 k8]; kl = [0 k2 k3 k5 k7 k9];
F_out = rN2in; R_gas_in = rN2in; Hi = [0 H.CO2 H.CH4 H.H2S H.NH3 H.H2O];
for i = 2:ng
    R_transfer  = -kla*(((GT(i)*R*T/Vg)/Hi(i))-C(kk(i)))*V; % mol/d    
    if i == 3
        R_transfer  = -kla*(((GT(i)*R*T/Vg)/Hi(i))-CT(kk(i)))*V; % mol/d
    end
    dMTdt(kl(i)) = dMTdt(kl(i))-R_transfer;
    R_gas_in = [R_gas_in R_transfer]; F_out = F_out + R_transfer;
end
%---------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------
% DEFINE DERIVATIVES IN GAS PHASE

for i = 1:ng
    R_gas_out = R_gas_in(i)-F_out*GT(i)/sum(GT);
    dMTdt = [dMTdt;R_gas_out];
end
%---------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------
% ADDITION OF OTHER DERIVATIVES

% Volume change
dMTdt = [dMTdt;Pm(8)]; 

% Cumulative CO2 and CH4
rCO2_out = F_out*GT(2)/sum(GT); dMTdt = [dMTdt;rCO2_out]; 
rCH4_out = F_out*GT(3)/sum(GT); dMTdt = [dMTdt;rCH4_out];

% additional output for rate information
R = [Rp.max(1)*ftot(1)*CT(1)*S(1,2);-1*Rp.max(2)*ftot(2)*CT(11)*S(2,2)];
end