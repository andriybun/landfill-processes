function [out] = Store_C_D(~,~,flag)
% function that is called by ODEsolver after each timestep
% Saves calculated Derived concentrations and times (C_D, t_D)

% Initialization of vectors to store C_D & t_D
global C_D t_D timer_flag C_D_o t_D_o 

if strcmp(flag,'init') == 1
    C_D_o = []; t_D_o = []; 
end
    
if strcmp(flag,'done')==0 && strcmp(flag,'init')==0
    C_D_o = [C_D_o;C_D'];
    t_D_o = [t_D_o t_D];
    % When timer_flag == 1, integration is flagged to stop
    out = timer_flag;
end

% Save C_D & t_D in mat file at the end of the integration
if strcmp(flag,'done') == 1
    % Round times to integer values & remove duplicate times
    t_D_o = round(t_D_o); [t_D_o, id] = unique(t_D_o);
    % Remove duplicate C_D
    nC = size(C_D_o); C_D_ot = zeros(length(t_D_o),nC(2));
    for i = 1:nC(2); C_D_ot(:,i) = C_D_o(id,i); end
    C_D_o = C_D_ot;
end
end