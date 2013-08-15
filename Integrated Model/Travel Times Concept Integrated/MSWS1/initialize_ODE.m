function [IC, Comp, Pm, S, Rp, H] = initialize_ODE(modus,xx)

%---------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------
% READ CSV FILE INTO XX
if modus == 0
    
    a = fopen('Pmatrix.csv'); N = 13; M = 38; XX = [];
    for i=1:M
        X = textscan(a,'%s',N);
        XX = [XX;X{:}'];
    end
    fclose('all');
    %---------------------------------------------------------------------------------------------------
    
    %---------------------------------------------------------------------------------------------------
    % DEFINE INITIAL CONCENTRATIONS AND PROCESS PARAMETERS
    
    % Initial concentrations of master species ((C-)mol/L)
    nM = 12; nc = 3; nini = 2; noutd = 5;
    
    Comp.master = XX(1,2:2+nM-1);
    Comp.masteri = str2double(XX(2,2:2+nM-1));
    
    Comp.const = XX(5,2:2+nc-1);
    Comp.consti = str2double(XX(6,2:2+nc-1));
    
    Comp.ini = XX(5,7:7+nini-1);
    Comp.inii = str2double(XX(6,7:7+nini-1));
    
    Comp.out =             [XX(33,2:end) XX(35,2:end) XX(37,2:N-noutd)];
    Comp.outi = str2double([XX(34,2:end) XX(36,2:end) XX(38,2:N-noutd)]);
    
    Comp.all = [Comp.master Comp.const Comp.ini Comp.out];
    Comp.alli = [Comp.masteri Comp.consti Comp.inii Comp.outi];
    
    % Stoichiomatrix
    nr = 7;
    S = str2double(XX(14:14+nr-1,2:2+nM-1));
    
    % Process parameters
    nP = 8;
    Pm = str2double(XX(10,2:2+nP-1));
    
    % Maximum rate constants (corrected for temperature)
    nrc = 1;
    Rp.max = str2double(XX(24:24+nr-1,2:2+nrc-1));
    
    T = Pm(1);     % K
    Rp.max(1) = Rp.max(1)*exp(64000/8.314*(1/303-1/T)); % Veeken & Hamelers 1998
    Rp.max(2) = Rp.max(2)*exp(64000/8.314*(1/303-1/T)); % Veeken & Hamelers 1998
    Rp.max(3) = Rp.max(3)*exp(64000/8.314*(1/303-1/T)); % Veeken & Hamelers 1998
    
    % Inhibition constants
    nri = 6;
    Rp.inhib = str2double(XX(24:24+nr-1,5:5+nri-1));
else
    Comp = xx{1}; Pm = xx{2}; S = xx{3}; Rp = xx{4};  T = Pm(1);    
end
%---------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------
%  DEFINE INITIAL H2CO3.TOT AND H+.TOT WITH INITIAL pCO2 AND pH

Ci = initialize_ORI(Comp,1);
k1 = find(strcmp('H2CO3.tot',Comp.all));
Comp.masteri(k1) = Ci(k1);  
k1 = find(strcmp('H+.tot',Comp.all));
Comp.masteri(k1) = Ci(k1);
Comp.alli = [Comp.masteri Comp.consti Comp.inii Comp.outi];
%---------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------
% SET INITIAL GAS CONCENTRATIONS ACCORDING TO HENRY'S LAW

R = Pm(3); p = Pm(2); Vg = Pm(5); 

H.CO2 = 1/0.027/0.9869*exp(-0.02629*(T-298)); % L*atm/mol 
H.CH4 = 1/0.0014/0.9869*exp(-0.01929*(T-298)); % L*atm/mol 
H.NH3 = 0.0164*exp(-0.0316*(T-298)); % L*atm/mol 
H.H2O = 0.00056835*exp(-0.0551*(T-298)); % L*atm/mol 
H.H2S = 10*exp(-0.0240*(T-298)); % L*atm/mol 

k1 = find(strcmp('H2CO3.con',Comp.all));
CO2i = Ci(k1)*H.CO2*Vg/R/T; % mol px = Hx*Cx nx = pxVg/RT
k1 = find(strcmp('NH3.con',Comp.all));
NH3i = Ci(k1)*H.NH3*Vg/R/T;
% NH3i = 0;
H2Oi = H.H2O*Vg/R/T; 
% H2Oi = 0; 
k1 = find(strcmp('H2S.con',Comp.all));
H2Si = Ci(k1)*H.H2S*Vg/R/T;
% H2Si =0;
k1 = find(strcmp('CH4',Comp.all));
CH4i = Ci(k1)*H.CH4*Vg/R/T;
N2i = p*Vg/R/T-CO2i-NH3i-H2Oi-H2Si-CH4i;
%---------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------
% STORE INITIAL CONCENTRATIONS AND PARAMETERS

Vli = Pm(4); MTi = Comp.masteri*Vli; GTi = [N2i CO2i CH4i H2Si NH3i H2Oi]; 
CO2_outi = 0; CH4_outi = 0; IC = [MTi GTi Vli CO2_outi CH4_outi];

if modus == 1
    if N2i < 0 
        IC = zeros(1,length(IC))+1e-6;
    end
end
end