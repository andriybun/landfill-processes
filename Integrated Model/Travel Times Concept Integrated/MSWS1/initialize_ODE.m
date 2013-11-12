function [IM, Comp, Pm, S, Rp] = initialize_ODE(fn)

%---------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------
% READ CSV FILE INTO XX

[MT_l, MT_g, OUTPUT, PARAMETERS, INHIBITION_TYPE] = read_Pmatrix(fn);

%---------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------
% DEFINE INITIAL CONCENTRATIONS AND PROCESS PARAMETERS

% Initial concentrations of master species ((C-)mol/L)
Comp.master = MT_l(1,5:end);
Comp.masteri = str2double(MT_l(2,5:end));

Comp.out = OUTPUT(2:end,1)';
Comp.outi = str2double(OUTPUT(2:end,2))';

Comp.all = [Comp.master Comp.out];
Comp.alli = [Comp.masteri Comp.outi];

% S
S = str2double(MT_l(3:end,5:end));

% Parameters
Pm.name = PARAMETERS(2:end,1);
Pm.value = str2double(PARAMETERS(2:end,2));

% Maximum rate constants (corrected for temperature)
Rp.max = str2double(MT_l(3:end,2:3));

T = Pm.value(1);     % K
k1 = find(strcmp('Yes', MT_l(3:end,4)));
Rp.max(k1,1) = Rp.max(k1,1).*exp(64000/8.314*(1/303-1/T)); % Veeken & Hamelers 1998

% Inhibition constants
Rp.inhib = str2double(INHIBITION_TYPE);

% Gas constants
Rp.gasv = str2double(MT_g(2:end,2:end));
Rp.gasn = MT_g(2:end,1);

%---------------------------------------------------------------------------------------------------

% %---------------------------------------------------------------------------------------------------
% %  DEFINE INITIAL TOTAL MASSES IN LIQUID AND GAS 
% 
k1 = strfind(Comp.master,'.logact'); k1 = find(cellfun(@isempty,k1)==0);
Mout_add = cellfun(@(x) [x(1:end-7) '.tot'], Comp.master(k1), 'UniformOutput', false);
k2 = strfind(Comp.master,'pH'); k2 = find(cellfun(@isempty,k2)==0);
H_add = cellfun(@(x) 'H+.tot', Comp.master(k2), 'UniformOutput', false);
Mout_add = [Mout_add H_add]; k1 = [k1 k2]; Comp.all = [Comp.all Mout_add]; Comp.alli = [Comp.alli 0.001*ones(1,length(Mout_add))];
Ci = initialize_ORI(Comp,1);
Comp.master(k1) = Mout_add; Comp.masteri(k1) = Ci(end+1-length(k1):end);
Comp.all = [Comp.master Comp.out]; Comp.alli = [Comp.masteri Comp.outi];

R = Pm.value(3); p = Pm.value(2); Vg = Pm.value(5); 
MT_g_i = Ci(Rp.gasv(2:end,3)).*(Rp.gasv(2:end,1)).*Vg./R./T; % mol px = Hx*Cx nx = pxVg/R
N2i = p*Vg/R/T-sum(MT_g_i);
MT_g_i = [N2i;MT_g_i];
%---------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------
% STORE INITIAL CONCENTRATIONS AND PARAMETERS
Vli = Pm.value(4); MT_l_i = Comp.masteri*Vli;

% Define cumulative outputs of total masses in liquid and gas
k1 = strfind(Comp.out,'.cum'); k1 = find(cellfun(@isempty,k1)==0);
Mcum_add = cellfun(@(x) x(1:end-4), Comp.out(k1), 'UniformOutput', false);
Rp.gascum = find(ismember(Rp.gasn,Mcum_add));  

IM = [MT_l_i MT_g_i' Vli Comp.outi(k1)];
end