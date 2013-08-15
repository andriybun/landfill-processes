%----------------------------------------------MSWS1------------------------------------------------    
function [simm] = integrate_bioreactor(P,ExtraPar)

% clc
% clear all
close all

pname = '../';
addpath(genpath(pname));
% modus = ExtraPar.modus;
modus = 0;
%-------------------------------Defining Parameters and Initial Conditions--------------------------

if modus == 0
    [IC, Comp, Pm, S, Rp, H] = initialize_ODE(modus,modus);
    ORI = initialize_ORI(Comp,0);
end

if modus == 1
    Pm = ExtraPar.Pm; Rp = ExtraPar.Rp; S = ExtraPar.S; IC = ExtraPar.IC; ORI = ExtraPar.ORI; 
    H = ExtraPar.H; tm = ExtraPar.tm; Comp = ExtraPar.Comp; T = Pm(1); Vli = Pm(4);
        
    % replacing parameters to be fitted by DREAM
%     % Case 1
%     IC(11) = P(1)*Vli;          IC(12) = P(2)*Vli;          Rp.inhib(2,1) = P(3);
%     Rp.max(1) = P(4)*exp(64000/8.314*(1/301-1/T));
%     Rp.max(2) = P(5)*exp(64000/8.314*(1/298-1/T));
%     Rp.max(3) = P(6)*exp(64000/8.314*(1/298-1/T));
%     Rp.inhib(1,3) = 10^-P(7);   Rp.inhib(2,3) = 10^-P(8);   Rp.inhib(1,4) = P(9);
%     Rp.inhib(2,5) = P(10);      Rp.inhib(2,6) = P(11);      Rp.inhib(3,6) = P(12);
%     Rp.max(7) = P(13);          Pm(6) = P(14);

%     % Case 2
%     IC(11) = P(1)*Vli;          IC(12) = P(2)*Vli;
%     Rp.max(7) = P(3);           Pm(6) = P(4);

      % case 4
      Comp.masteri(1) = P(1); Comp.masteri(2) = P(2); Comp.masteri(4) = P(3); Comp.masteri(5) = P(4); Comp.masteri(7) = -2*(P(11)-P(5))-(P(12)-P(13));
      Comp.masteri(8) = P(5); Comp.masteri(9) = P(6); Comp.masteri(10) = P(7); Comp.masteri(11) = P(8); Comp.masteri(12) = P(9);
      Comp.inii(1) = P(10); Comp.consti(1) = P(11); Comp.consti(2) = P(12); Comp.consti(3) = P(13);
      
      xx = {Comp, Pm, S, Rp};
      [IC, Comp, Pm, S, Rp, H] = initialize_ODE(modus,xx);
      ORI = initialize_ORI(Comp,0);
      
end
%------------------------------------------Simulation-----------------------------------------------

% global variables to store intermediate ORCHESTRA results
global Call V2 tt Rall
Call = [];  V2 = []; tt = []; Rall = [];

% Solve Differential Equation
if modus == 0
    options = odeset('OutputFcn','StoreC','Refine',1, 'AbsTol', 1e-8, 'RelTol', 1e-8);
end
if modus == 1
    options = odeset('OutputFcn','StoreC','Refine',1, 'AbsTol', 1e-4, 'RelTol', 1e-4);
end
trange = 0:1:800;
[t,MT] = ode15s(@bioreactor, trange, IC, options, Comp, Pm, S, Rp, ORI, H);

% safety measure for incomplete integration
if modus == 1 && length(t) < length(trange)
    t = trange;
    tt = trange;
    MT = zeros(length(trange),length(IC));
    Call = zeros(length(trange),length(Comp.outi))';
end

% ---------------------------------------Data-processing--------------------------------------------

[sim, ts] = plot_integration(modus, Comp, Pm, S, t, MT, Call, V2, tt, Rall);

if modus == 1
    % Sort simulated data
    simm = (sort_simdata(sim,ts,tm))';
end
end