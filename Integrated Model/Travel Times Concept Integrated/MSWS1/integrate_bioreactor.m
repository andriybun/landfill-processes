function [out, sim] = integrate_bioreactor(P,Extra)

% close all;
pname = '../MSWS1'; addpath(genpath(pname));

%% Run MSWS1 

switch Extra.sw
    case 0  % single run 

        tic
        [IM, Pm, ORI, Meas, ~] = initialize(Extra.fn);
        toc
        
        tic
        global C_D_o t_D_o
        options = odeset('NonNegative',[1:6 8:23],'OutputFcn','Store_C_D','Refine',1, 'AbsTol', 1e-5, 'RelTol',1e-5);
        trange = Pm.tstart:Pm.tstep:Pm.tend;
        [t,MT] = ode15s(@bioreactor, trange, IM, options, Pm, ORI);
        toc
        
        tic
        % Dataprocessing
        plot_integration(t, MT, t_D_o, C_D_o, Meas, Pm);
        out = 0; sim = 0;
        toc
        
    case 1  % DREAM run
        
        tic
        Pm = Extra.Pm; 
        % Replaces parameters with DREAM parameters
        Pm.Kmax(1) = P(1); Pm.Kmax(2) = P(2); Pm.Kmax(3) = P(3); Pm.Kmax(7) = P(4);
        Pm.Ki(2) = P(5); Pm.Ki(3) = P(6); Extra.IM(11) = P(7)*Pm.V_l; Extra.IM(12) = P(8)*Pm.V_l;
        % Also correct decay parameters!
        Pm.Kmax(5,1) = 0.05*P(2);%*0.15; 
        Pm.Kmax(6,1) = 0.05*P(3);%*0.094;
        % Get sigma parameters
        Pm.sigma = P(9:end)';
        
        global C_D_o t_D_o
        options = odeset('NonNegative',[1:6 8:23],'OutputFcn','Store_C_D','Refine',1, 'AbsTol', 1e-5, 'RelTol',1e-5);
        trange = Pm.tstart:Pm.tstep:Pm.tend;
        [t,MT] = ode15s(@bioreactor, trange, Extra.IM, options, Pm, Extra.ORI);
        % Safety measure for incomplete integration & negative concentrations VFA, NH4
        % Time is set as initial range and all concentrations are set to zero
        if length(t) < length(trange) || min(min([MT(:,2) MT(:,4)])) < 0
            t = trange; t_D_o = trange;
            MT = zeros(length(trange),length(Extra.IM)); C_D_o = zeros(length(trange),35);
        end
                        
        % Dataprocessing
%         plot_integration(t, MT, t_D_o, C_D_o, Extra.Meas, Pm);
        sim.t = t; sim.MT = MT; sim.t_D_o = t_D_o; sim.C_D_o = C_D_o;
        out = post_process(t, MT, t_D_o, C_D_o, Extra.Meas, Pm);
        toc
end
end