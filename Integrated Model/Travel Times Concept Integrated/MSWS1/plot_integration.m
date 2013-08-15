function [sim, ts] = plot_integration(modus, Comp, Pm, S, t, MT, Call, V2, tt, Rall)

R = Pm(3); T = Pm(1); p = Pm(2); Vg = Pm(5);
%---------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------
% DEFINE EXPERIMENTAL DATA

if modus == 0
    Meas = load_data(Pm);
end
%---------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------
% RE-DEFINE SIMULATED DATA

[nreac,ncomp] = size(S);    ng = 6;
GT = MT(:,ncomp+1:ncomp+ng); V = MT(:,ncomp+ng+1); CO2_out = MT(:,ncomp+ng+2); CH4_out = MT(:,ncomp+ng+3);
MT = MT(:,1:ncomp);

Bio_s = (CO2_out+CH4_out)/1000;

k1 = find(strcmp('H[Acetate].tot',Comp.master));
VFA_s = MT(:,k1)./V;

k1 = find(strcmp('H+.con',Comp.out));
pH_s = -log10(Call(:,k1));

k1 = find(strcmp('NH3.tot',Comp.master));
NH4_s = MT(:,k1)./V;

CO2p_s = GT(:,2)*R*T/Vg;
CH4p_s = GT(:,3)*R*T/Vg;

sim = {Bio_s pH_s' VFA_s NH4_s CH4p_s CO2p_s};
ts = {t tt t t t t};

% Calculate Ionic strength
I = 0.5.*(sum(Call(:,8:19),2)+4*sum(Call(:,20:22),2));

% Calculate EC from individual concentrations
EC_i = [349.8 38.6 119/2/2 76.3 73.5 65 44.5 44.5/2 198.6 80/2/2 52 50.1 160 119 138.6];
C_i = Call(:,8:22); EC_i = C_i*EC_i'./1000;

% Calculate EC from correlation with Ionic strength considering ECm, pHm, OHm, NH4m, VFAm
EC_s = 40.6870.*I+6.4450; 

% % Calculate EN from all concentrations
% EN = Call;
% charge = [0 0 0 0 0 0 0 1 -1 1 -1 1 -1 -1 -1 -1 -1 -1 1 -2 2 -2 1 2 1 2 -2 0 0 0 -1];
% EN = EN*charge';
% 
% H_tot = []; SO4_tot = [];
% for i = 1:length(tt)
%     k1 = find(tt(i) == t);
%     H_tot(i) = MT(k1,7)/Pm(4);
%     SO4_tot(i) = MT(k1,8)/Pm(4);
% end
% EN = EN - (H_tot'+ 2.1 + (0.14-2*SO4_tot'));
%---------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------
% FIGURE 1: MASTER SPECIES IN TIME

if modus == 0
    figure(1)
    for i = 1:ncomp,
        subplot(3,4,i)
        hold on
        plot(t,MT(:,i)./V,'-r');
        title (Comp.master{i});
    end
end
%---------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------
% FIGURE 2: CONCENTRATIONS IN TIME

if modus == 0 
    figure(2)
    for i = 1:length(Call(1,:)),
        subplot(6,6,i)
        hold on
        plot(tt,Call(:,i),'-b');
        title (Comp.out{i});
    end
end
%---------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------
% FIGURE 3: SIMULATED AND EXPERIMENTAL BIOGAS, VFA AND pH IN TIME

if modus == 0 
    h = figure(3);
    set(h,'Position',[188 388 1303 300]);
    subplot(1,3,1)
    hold on
    plot(t,Bio_s*R*T/p,'r-','LineWidth',5)
    plot(Meas.Bio_tm,Meas.Bio_m*R*T/p,'x','MarkerSize',7)
    xlabel('Time (days)','FontSize',15)
    ylabel('CO2+CH4 (m^{3})','FontSize',15)
    title('Biogas','FontSize',15);
    set(gca,'FontSize',15)
    
    subplot(1,3,2)
    plot(t,VFA_s,'r-','LineWidth',5)
    hold on
    plot(Meas.VFA_tm,Meas.VFA_m,'x','MarkerSize',7)
    xlabel('Time (days)','FontSize',15)
    ylabel('VFA (mol/L)','FontSize',15)
    title('Total VFA','FontSize',15);
    set(gca,'FontSize',15)

    subplot(1,3,3)
    plot(tt,pH_s,'r-','LineWidth',5)
    hold on
    plot(Meas.pH_tm,Meas.pH_m,'x','MarkerSize',7)
    xlabel('Time (days)','FontSize',15)
    ylabel('pH','FontSize',15)
    title('pH','FontSize',15);
    set(gca,'FontSize',15)
end
%---------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------
% FIGURE 4: SIMULATED AND EXPERIMENTAL NH4, EC, pCO2 AND pCH4 IN TIME

if modus == 0
    figure(4)
    subplot(2,2,1)
    plot(Meas.NH4_tm,Meas.NH4_m,'xb')
    hold on
    plot(t,MT(:,4)./V,'-xr')
    title('NH4')
    xlabel('Time (days)')
    ylabel('NH4')
    
    subplot(2,2,2)
    plot(Meas.EC_tm,Meas.EC_m,'xb')
    hold on
    plot(tt,EC_s,'-xr')
    plot(tt,EC_i,'-xm')
    title('EC')
    xlabel('Time (days)')
    ylabel('EC')
    
    subplot(2,2,3)
    plot(Meas.CH4p_tm,Meas.CH4p_m,'xb')
    hold on
    plot(t,CH4p_s,'-xr')
    title('pCH4')
    xlabel('Time (days)')
    ylabel('pCH4')
    
    subplot(2,2,4)
    plot(Meas.CO2p_tm,Meas.CO2p_m,'xb')
    hold on
    plot(t,CO2p_s,'-xr')
    title('pCO2')
    xlabel('Time (days)')
    ylabel('pCO2')
end
%---------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------
% FIGURE 5: SIMULATED AND EXPERIMENTAL DATA FOR PRESENTATION

if modus == 0
    h = figure(5);
    set(h,'Position',[10 10 1000 700]);
    subplot(2,3,1)
    hold on
    plot(t,Bio_s*R*T/p,'r-','LineWidth',5)
    plot(Meas.Bio_tm,Meas.Bio_m*R*T/p,'x','MarkerSize',7)
    xlabel('t (d)','FontSize',15)
    ylabel('CO_2&CH_4 (m^{3})','FontSize',15)
    title('Biogas','FontSize',15);
    set(gca,'FontSize',15)
    
    subplot(2,3,2)
    plot(t,VFA_s,'r-','LineWidth',5)
    hold on
    plot(Meas.VFA_tm,Meas.VFA_m,'x','MarkerSize',7)
    xlabel('t (d)','FontSize',15)
    ylabel('VFAx (M)','FontSize',15)
    title('Total VFAx','FontSize',15);
    set(gca,'FontSize',15)

    subplot(2,3,3)
    plot(tt,pH_s,'r-','LineWidth',5)
    hold on
    plot(Meas.pH_tm,Meas.pH_m,'x','MarkerSize',7)
    xlabel('t (d)','FontSize',15)
    ylabel('pH','FontSize',15)
    title('pH','FontSize',15);
    set(gca,'FontSize',15)
    
    subplot(2,3,4)
    plot(Meas.NH4_tm,Meas.NH4_m,'xb')
    hold on
    plot(t,MT(:,4)./V,'-xr')
    title('NH_4^+&NH_3 (M)','FontSize',15)
    xlabel('t (d)','FontSize',15)
    ylabel('NH_4_&NH_3 (M)','FontSize',15)
    set(gca,'FontSize',15)
    
    subplot(2,3,5)
    plot(Meas.CH4p_tm,Meas.CH4p_m,'xb')
    hold on
    plot(t,CH4p_s,'-xr')
    title('pCH_4','FontSize',15)
    xlabel('t (d)','FontSize',15)
    ylabel('pCH_4 (atm)','FontSize',15)
    set(gca,'FontSize',15)
    
    subplot(2,3,6)
    plot(Meas.CO2p_tm,Meas.CO2p_m,'xb')
    hold on
    plot(t,CO2p_s,'-xr')
    title('pCO_2','FontSize',15)
    xlabel('t (d)','FontSize',15)
    ylabel('pCO_2 (atm)','FontSize',15)
    set(gca,'FontSize',15)
end

%---------------------------------------------------------------------------------------------------

% %---------------------------------------------------------------------------------------------------
% % FIGURE 6: HISTOGRAM P_fit_wide
% 
% if modus == 0
%     
% %     load Optimization_P_fit_wide
%     load Optimization
% %     Pnames = {'Ci_x(meth) (M)' 'Ci_x(sulf) (M)' 'k(NH_3) (d^{-1})' 'k_l_a (d^{-1})'};
%     Pnames = {'Ci_{x,meth} (mM)' 'Ci_{x,sulf} (M)'  'K_{s,meth} (M)' 'k_{hyd} (d^{-1})' 'q_{s,meth}^{max}(d^{-1})' ...
%                'q_{s,sulf}^{max}(d^{-1})' 'pH_{i,hyd}' 'pH_{i,meth}' 'K_{i,VFAx} (M)' 'K_{i,NH3} (M)' ...
%                'K_{i,H2S,meth} (M)'   'K_{i,H2S,sulf} (M)'  'rNH_3 (d^{-1}x10^{-3})' 'k_l_a (d^{-1})'};
%     Zind=max(find(Z(:,1)>0)); %the end-value of Z
%     [minZ minZidx] = max(Z(MCMCPar.m0+1:Zind,MCMCPar.n+1));
%     minZidx=minZidx+MCMCPar.m0;  %the index corresponding to the parameter set with the highest probability
%     
%     figure(6)
%     for i = 1:MCMCPar.n
%         subplot(4,4,i);
%         if i == 1 || i == 13
%             Z(round(0.8.*Zind):Zind,i) = Z(round(0.8.*Zind):Zind,i)*1000;
%         end
%         hist(Z(round(0.8.*Zind):Zind,i),10)
%         
%         xlabel(Pnames{i},'FontSize',14)
%         ylabel('Frequency','FontSize',12)
%         set(gca,'FontSize',14);
% %             xlim([ParRange.minn(i) ParRange.maxn(i)])
%         h = findobj(gca,'Type','patch');
%         set(h,'FaceColor',[0.5 0.5 0.5])
%     end
% end

%---------------------------------------------------------------------------------------------------

% %---------------------------------------------------------------------------------------------------
% % FIGURE 7: BALANCES
% 
% if modus == 0
%     
%     figure(7)
%     plot(tt,EN);
%     title('EN Balance','FontSize',15)
%     xlabel('t (d)','FontSize',15)
%     ylabel('EN (charge)','FontSize',15)
%     set(gca,'FontSize',15)
% end
%---------------------------------------------------------------------------------------------------

% %---------------------------------------------------------------------------------------------------
% % FIGURE 8: RATES
% 
% if modus == 0
%     
%     figure(8)
%     plot(tt,Rall(:,1),'b')
%     hold on
%     plot(tt,Rall(:,2),'r')
%     title('Rates','FontSize',15)
%     xlabel('t (d)','FontSize',15)
%     ylabel('Rates (mol d^{-1})','FontSize',15)
%     legend('R(hyd)' ,'R(meth)')
% end

%---------------------------------------------------------------------------------------------------
% FIGURE 9: GAS IN TIME

GT_names = {'N2' 'CO2' 'CH4' 'H2S' 'NH3' 'H2O'};
if modus == 0
    figure(9)
    for i = 1:ng,
        subplot(3,4,i)
        hold on
        plot(t,GT(:,i)./Vg,'-r');
        title (GT_names{i});
    end
end
%---------------------------------------------------------------------------------------------------
end