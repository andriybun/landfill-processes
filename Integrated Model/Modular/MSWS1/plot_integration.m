function plot_integration(t, MT, t_D, C_D, Meas, Pm)

% Define total concentrations in liquid phase (MT_l) and gas phase (MT_g) 
MT_l = MT(:,1:Pm.ncompl); MT_g = MT(:,Pm.ncompl+1:Pm.ncompl+Pm.ncompg);

% FIGURE 1: Plots all Total concentrations in liquid and gas phase and cumulative outputs
% Calculates proper number of subplots
a = floor(sqrt(length(MT(1,:)))); b = a;
if a*b < length(MT) || a*b == 1
    b = b+1;
    if a*b < length(MT(1,:))
        a = a+1;
    end
end
figure;
CT = [MT_l./Pm.V_l MT_g MT(:,end-1:end)];
for i = 1:length(MT(1,:)),
    subplot(a,b,i)
    hold on
    plot(t,CT(:,i),'-r');
    title (Pm.names_CT{i});
end

% FIGURE 2: Plots all Derived concentrations
% Calculates proper number of subplots
a = floor(sqrt(length(C_D(1,:)))); b = a;
if a*b < length(C_D(1,:)) || a*b == 1
    b = b+1;
    if a*b < length(C_D(1,:))
        a = a+1;
    end
end
figure;
for i = 1:length(C_D(1,:)),
    subplot(a,b,i)
    hold on
    plot(t_D,C_D(:,i),'-b');
    title (Pm.names_CD{i});
end

% Dataspecific figures
if Pm.dataset > 0
    
    % Names of modelled concentrations, transformed into the proper units
    CO2_out = MT(:,Pm.ncompl+Pm.ncompg+1);
    CH4_out = MT(:,Pm.ncompl+Pm.ncompg+2);
    Biogas_sim = (CO2_out+CH4_out)*Pm.R_gas*Pm.T/Pm.p/1000;  % m^3
    VFA_sim = MT_l(:,2)./Pm.V_l; % mol/L
    pH_sim = -log10(C_D(:,25));
    NH4_sim = MT_l(:,4)./Pm.V_l; % mol/L
    CO2p_sim = MT_g(:,2)*Pm.R_gas*Pm.T/Pm.V_g; % atm
    CH4p_sim = MT_g(:,3)*Pm.R_gas*Pm.T/Pm.V_g; % atm
    
    % Plots cumulative biogas (simulated & experiment)
    h = figure;
    set(h,'Position',[10 10 1000 700]);
    subplot(2,3,1)
    hold on
    plot(t,Biogas_sim,'xr-','LineWidth',5)
    plot(Meas.Biogas_tm,Meas.Biogas_m,'x','MarkerSize',7)
    xlabel('t (d)','FontSize',15)
    ylabel('CO_2&CH_4 (m^{3})','FontSize',15)
    title('Biogas','FontSize',15);
    set(gca,'FontSize',15)
    
    % Plots VFA (simulated & experiment)
    subplot(2,3,2)
    plot(t,VFA_sim,'xr-','LineWidth',5)
    hold on
    plot(Meas.VFA_tm,Meas.VFA_m,'x','MarkerSize',7)
    xlabel('t (d)','FontSize',15)
    ylabel('VFAx (M)','FontSize',15)
    title('Total VFAx','FontSize',15);
    set(gca,'FontSize',15)
    
    % Plots pH (simulated & experiment)
    subplot(2,3,3)
    plot(t_D,pH_sim,'xr-','LineWidth',5)
    hold on
    plot(Meas.pH_tm,Meas.pH_m,'x','MarkerSize',7)
    xlabel('t (d)','FontSize',15)
    ylabel('pH','FontSize',15)
    title('pH','FontSize',15);
    set(gca,'FontSize',15)
    
    % Plots NH4 + NH3 (simulated & experiment)
    subplot(2,3,4)
    plot(Meas.NH4_tm,Meas.NH4_m,'xb')
    hold on
    plot(t,NH4_sim,'-xr')
    title('NH_4^+&NH_3 (M)','FontSize',15)
    xlabel('t (d)','FontSize',15)
    ylabel('NH_4_&NH_3 (M)','FontSize',15)
    set(gca,'FontSize',15)
    
    % Plots pCH4 (simulated & experiment)
    subplot(2,3,5)
    plot(Meas.CH4p_tm,Meas.CH4p_m,'xb')
    hold on
    plot(t,CH4p_sim,'-xr')
    title('pCH_4','FontSize',15)
    xlabel('t (d)','FontSize',15)
    ylabel('pCH_4 (atm)','FontSize',15)
    set(gca,'FontSize',15)
    
    % Plots pCO2 (simulated & experiment)
    subplot(2,3,6)
    plot(Meas.CO2p_tm,Meas.CO2p_m,'xb')
    hold on
    plot(t,CO2p_sim,'-xr')
    title('pCO_2','FontSize',15)
    xlabel('t (d)','FontSize',15)
    ylabel('pCO_2 (atm)','FontSize',15)
    set(gca,'FontSize',15)
end
end