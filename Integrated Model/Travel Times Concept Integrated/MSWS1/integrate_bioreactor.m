function [] = integrate_bioreactor1(fn)

close all; pname = '../../Model'; addpath(genpath(pname));

%% Set Dataset and modus
modus = 0;

%% Run MSWS1 

[IM, Comp, Pm, S, Rp] = initialize_ODE(fn);
ORI = initialize_ORI(Comp,0);

% P = [ 0.0724   0.581      3.368       0.0058      0.0001      0.3484    0.0401    0.0097];
% % replacing parameters to be fitted by DREAM
% Rp.max(1,1) = P(1); Rp.max(2,1) = P(2); Rp.max(3,1) = P(3); Rp.max(7,1) = P(4);
% Rp.inhib(2,4) = P(5); Rp.inhib(3,4) = P(6);
% Comp.masteri(11) = P(7); IM(11) = P(7)*IM(22);
% Comp.masteri(12) = P(8); IM(12) = P(8)*IM(22);
% % Also correct decay parameters!
% Rp.max(5,1) = 0.05*P(2)*S(2,11);    Rp.max(6,1) = 0.05*P(3)*S(3,12);

options = odeset('OutputFcn','Store_Orchestra_Results','Refine',1, 'AbsTol', 1e-8, 'RelTol', 1e-8);
trange = Pm.value(strcmp(Pm.name,'tstart')):Pm.value(strcmp(Pm.name,'tstep')):Pm.value(strcmp(Pm.name,'tend'));
tic
[t,MT] = ode15s(@bioreactor, trange, IM, options, Comp, Pm, S, Rp, ORI);
toc

% Dataprocessing
load Results_Orchestra; V_o = MS_o(:,end)*ones(1,length(Comp.outi)); MS_o = MS_o(:,1:end-1); CS_o = MS_o./V_o;
[simr, ts] = plot_integration(Comp, Pm, S, t, MT, CS_o, tS_o, Rp, modus);
        
end