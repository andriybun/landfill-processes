function out = post_process(t, MT, t_D, C_D, Meas, Pm)

% Define simulated data
CO2_out = MT(:,Pm.ncompl+Pm.ncompg+1);
CH4_out = MT(:,Pm.ncompl+Pm.ncompg+2);
Biogas_sim = (CO2_out+CH4_out)*Pm.R_gas*Pm.T/Pm.p/1000;  % m^3
VFA_sim = MT(:,2)./Pm.V_l; % mol/L
NH4_sim = MT(:,4)./Pm.V_l; % mol/L
CO2p_sim = MT(:,17)*Pm.R_gas*Pm.T/Pm.V_g; % atm
CH4p_sim = MT(:,18)*Pm.R_gas*Pm.T/Pm.V_g; % atm

% Derived concentration data is interpolated to get data for each day
if max(C_D(:,25)) == 0; pH_sim = interp(C_D(:,25),150); else pH_sim = interp(-log10(C_D(:,25)),150); end 
t_D_pH = interp(t_D,150)'; [t_D_pH, id] = unique(round(t_D_pH)); pH_sim = pH_sim(id);

% Sort simulated data based on experimental times & define sigma
id = ismember(t,Meas.Biogas_tm); Biogas_sim = Biogas_sim(id); 
id = ismember(t,Meas.VFA_tm); VFA_sim = VFA_sim(id); 
id = ismember(t_D_pH,Meas.pH_tm); pH_sim = pH_sim(id); 
id = ismember(t,Meas.NH4_tm); NH4_sim = NH4_sim(id); 
id = ismember(t,Meas.CO2p_tm); CO2p_sim = CO2p_sim(id); 
id = ismember(t,Meas.CH4p_tm); CH4p_sim = CH4p_sim(id); 

% Define sigma 
if length(Pm.sigma) == 1
    set_sigma = Pm.sigma;
elseif length(Pm.sigma) == 6
    s_Biogas(1:length(Biogas_sim)) = Pm.sigma(1); s_VFA(1:length(VFA_sim)) = Pm.sigma(2);
    s_pH(1:length(pH_sim)) = Pm.sigma(3); s_NH4(1:length(NH4_sim)) = Pm.sigma(4);
    s_CO2p(1:length(CO2p_sim)) = Pm.sigma(5); s_CH4p(1:length(CH4p_sim)) = Pm.sigma(6);
    set_sigma = [s_Biogas s_VFA s_pH s_NH4 s_CO2p s_CH4p]';
end

% Calculate ln_p based on multivariate gaussian distribution of error --> e = di - Fi
% set_exp = [Meas.Biogas_m; Meas.VFA_m; Meas.pH_m; Meas.NH4_m; Meas.CO2p_m; Meas.CH4p_m];
% set_sim = [Biogas_sim; VFA_sim; pH_sim; NH4_sim; CO2p_sim; CH4p_sim];
% try to put extra weight on biogas & pH
W = 10;
set_exp = [Meas.Biogas_m.*W; Meas.VFA_m; Meas.pH_m.*W; Meas.NH4_m; Meas.CO2p_m; Meas.CH4p_m];
set_sim = [Biogas_sim.*W; VFA_sim; pH_sim.*W; NH4_sim; CO2p_sim; CH4p_sim];

N = [length(Meas.Biogas_m); length(Meas.VFA_m); length(Meas.pH_m); length(Meas.NH4_m); length(Meas.CO2p_m); length(Meas.CH4p_m)];

% p = (2*pi()*set_sigma^2)^(-0.5*N)*exp(-0.5*set_sigma^-2*sum((set_exp-set_sim).^2));
ln_p = -0.5*sum(N)*log(2*pi()) + sum(-0.5.*N.*log(Pm.sigma.^2)) + sum(-0.5.*set_sigma.^-2.*(set_exp-set_sim).^2);
out = ln_p;
end



