function dMTdt = bioreactor(t, MT, Pm, ORI)

% Define total concentrations in liquid phase (MT_l) and gas phase (MT_g) 
MT_l = MT(1:Pm.ncompl); MT_g = MT(Pm.ncompl+1:Pm.ncompl+Pm.ncompg);

% bacterial concentrations below a limit are not active anymore for numerical stability
MT_l_b = MT_l(end-2:end); MT_l_b(MT_l_b < 1e-7) = 0; MT_l(end-2:end) = MT_l_b; 

%% Calculation of Derived concentrations with ORCHESTRA
global C_D t_D
% Replace negative values of MT_l
MT_o = MT_l; id = MT_o < 1e-7; id(Pm.Htoti) = 0; MT_o(id) = 1e-7;
% Calculate derived concentrations C_D with times t_D
CT = MT_o./Pm.V_l; ro = ORI.Calculate(1:Pm.ncompl, CT);
C_D = ro(Pm.ncompl+1:end); t_D = t; 
% Join total concentrations liquid phase and derived concentrations
C_all = [MT_l./Pm.V_l;C_D];
% Check for proper ORCHESTRA calculation. If Nan is present, calculation not valid --> integrations stopped. 
if max(isnan(ro)) > 0; dMTdt = zeros(length(MT),1); return; end;

%% Inhibition
% Calculate total inhibition factors I
I = inhibition(Pm, C_all);

%% Total Biochemical rates
% Calculate total biochemical rates Rbio, C_Ri specifies rate determining concentrations
Rbio = (Pm.Kmax.*I.*MT_l(Pm.C_Ri))'*Pm.S;

%% Mass transfer liquid/gas
% Calculate total mass transfer for liquid phase (Ftransfer_l) and gas phase (Ftransfer_g)
% Ftransfer_g == Ftransfer_l without zero entries
Ftransfer_l = zeros(Pm.ncompl,1);
Ftransfer_g  = -Pm.kla.*((MT_g.*Pm.R_gas.*Pm.T./Pm.V_g./Pm.H)-C_all(Pm.Ci_tf)).*Pm.V_l; % mol/d
Ftransfer_l(Pm.C_Ti_tf) = Ftransfer_g;

%% Mass transport gas
% Calculate total mass transport out of the gas phase
Ftransport = MT_g./sum(MT_g).*sum(Ftransfer_g);

%% Change of total concentrations liquid/gas
% check time of integration for very slow 'stuck' integrations, threshold is max time allowed (s) 
threshold = 20; sw = check_duration(t,threshold);
switch sw
    case 0
        % Definitition of total changes in total concentrations in liquid and gas phase. Additionally, cumulative outputs are defined.  
        dMTdt = [Rbio' - Ftransfer_l; ...   % liquid phase
            Ftransfer_g - Ftransport; ... % gas phase
            Ftransport(Pm.cumgasi)];    % cumulative gas production
    case 1
        % In case of 'stuck' or wrong integration, fast ending.
        dMTdt = zeros(length(MT),1);
end
end 

function I = inhibition(Pm, C_all)

% Define inhibiting concentrations
Cinh = C_all(Pm.Cinhi);
% Calculation of all inhibition factors per inhibtion mechanism per reaction (I_i)
I_i(Pm.finh==1) = Cinh(Pm.finh==1)./(Cinh(Pm.finh==1) + Pm.Ki(Pm.finh==1));
I_i(Pm.finh==2) = Pm.Ki(Pm.finh==2)./(Cinh(Pm.finh==2) + Pm.Ki(Pm.finh==2));
I_i(Pm.finh==3) = Pm.Ki(Pm.finh==3).^2./(Cinh(Pm.finh==3).^2 + Pm.Ki(Pm.finh==3).^2);
I_i(Pm.finh==4) = (1-Cinh(Pm.finh==4)./Pm.Ki(Pm.finh==4)).^0.401;
I_i(Pm.finh==5) = (1+(Cinh(Pm.finh==5)./Pm.Ki(Pm.finh==5)).^1.08).^-1;

% Deletion of imaginary numbers, mechanisms f4 & f5 might give these at certain inhibitor concentrations. 
% In this case their factor is 1. 
Pm.Rinh = Pm.Rinh(I_i == real(I_i));
I_i = I_i(I_i == real(I_i));

% All inhibition factors acting on the same reaction are multiplied per reaction 
I = ones(length(Pm.Kmax),1);
I(1) = prod(I_i(Pm.Rinh == 1));
I(2) = prod(I_i(Pm.Rinh == 2));
I(3) = prod(I_i(Pm.Rinh == 3));
end

function [sw] = check_duration(t, threshold)
%% Checks the duration of the integration. 
% When integration gets 'stuck', integration is forced to stop after threshold.
% timer_flag is a flag used in the function Store_C_D called by the ODE_solver after each step.
% If timer_flag is 1, integration is stopped, if sw is 1 all dMTdt are 0.
global timer_flag; persistent timer_0;

if t == 0
    timer_0 = tic;
end

tpass = toc(timer_0);
if tpass > threshold
    timer_flag = 1; sw = 1;
else
    timer_flag = 0; sw = 0;
end
end
