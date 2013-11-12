function dMTdt = bioreactor(t, MT, Comp, Pm, S, Rp, ORI)
%% Some vector definitions 

[~,ncompl] = size(S); ncompg = length(Rp.gasv(:,1)); MT_l = MT(1:ncompl); MT_g = MT(ncompl+1:ncompl+ncompg); V = MT(ncompl+ncompg+1);

D = MT_l(Rp.max(:,2)); K = Rp.max(:,1); L = ones(1,ncompl);

T = Pm.value(1); p = Pm.value(2); R_gas = Pm.value(3); Vg = Pm.value(5); rN2in = Pm.value(6)*1e-3*p/R_gas/T;

global MS tS

%% Measures to stabilize the calculations

    %% Make MT without negative concentrations for ORCHESTRA
    MT_o = MT_l;
    id = find(MT_o < 1e-10);
    id(id==find(strcmp('H+.tot',Comp.master))) = [];
    MT_o(id) = 1e-10;

    %% Force Bacterial masses to stay at 0 when very low
    id = find(MT_l < 1e-10);
    id(id < 10) = [];
    MT_l(id) = 0;

    %% Checks the duration of the integration. When integration gets 'stuck', integration is forced to stop after threshold. 
    global timer_0 timer_flag
    threshold = 60; % s
    if t > 0.001
        tpass = toc(timer_0);
        if tpass > threshold;
           timer_flag = 1;
        else
            timer_flag = 0;
        end
    else
        timer_flag = 0;
    end

%% Define vector of derivatives

switch timer_flag
    case 0
        %% Calculation of specific mass/concentrations with ORCHESTRA
        n = (1:ncompl); CT = MT_o./V; 
        ro = ORI.Calculate(n, CT);
        CS = ro(ncompl+1:end);
        MS = [CS.*V;V]; tS = t;
        
        %% Inhibition
        I = inhibition(Rp, Comp, MT_l, CS, V);

        %% (Bio)chemical reaction part of derivatives in liquid phase R
        R = sum((S.*(K.*I.*D*L)),1)';
        dMTdt = R; 

        %% Update derivatives in liquid phase for mass transfer H
        CS = [(MT_l./V);CS]; H_p = Rp.gasv(2:end,1); kla = Rp.gasv(2:end,2);
        H  = kla.*((MT_g(2:end).*R_gas.*T./Vg./H_p)-CS(Rp.gasv(2:end,3))).*V; % mol/d
        dMTdt(Rp.gasv(2:end,4)) = dMTdt(Rp.gasv(2:end,4)) + H;
                     
        %% Define derivatives for gas phase H - F  
        H = [rN2in; -H]; 
        F = MT_g./sum(MT_g).*sum(H);
        dMTdt = [dMTdt;H-F];
                    
        %% Volume change
        dMTdt = [dMTdt;Pm.value(7)];
        
        %% Cumulative CO2 and CH4
        dMTdt = [dMTdt;F(Rp.gascum)];
        
    case 1
        %% Define zero derivatives when t > threshold time
        dMTdt = zeros(length(MT),1);
end 
end