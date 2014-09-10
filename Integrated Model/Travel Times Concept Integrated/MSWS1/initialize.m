function [IM, Pm, OrchestraObj, Meas, Comp] = initialize(fn)

% Read CSV file into XX
[MT_l, MT_g, OUTPUT, PARAMETERS, INHIBITION_TYPE] = read_Pmatrix(fn);

% Initial concentrations of master species for ORCHESTRA ((C-)mol/L)
Comp.master = MT_l(1,5:end);
Comp.masteri = str2double(MT_l(2,5:end));
Comp.out = OUTPUT(2:end,1)';
Comp.outi = str2double(OUTPUT(2:end,2))';
Comp.all = [Comp.master Comp.out];
Comp.alli = [Comp.masteri Comp.outi];

% Stoichiometry matrix S
Pm.S = str2double(MT_l(3:end,5:end));

% Inhibition constants
Pm.Ki = str2double(INHIBITION_TYPE(2:end,4));
Pm.finh = str2double(INHIBITION_TYPE(2:end,1));
Pm.Rinh = str2double(INHIBITION_TYPE(2:end,2));
Pm.Cinhi = str2double(INHIBITION_TYPE(2:end,3));

% Gas constants
Pm.H = str2double(MT_g(2:end,2));
Pm.kla = str2double(MT_g(2:end,3));
Pm.Ci_tf = str2double(MT_g(2:end,4));
Pm.C_Ti_tf = str2double(MT_g(2:end,5));
Pm.gas_name = MT_g(2:end,1);

% Process parameters & others
Pm.T = str2double(PARAMETERS(strcmp('T',PARAMETERS(:,1)),2));
Pm.p = str2double(PARAMETERS(strcmp('p',PARAMETERS(:,1)),2));
Pm.R_gas = str2double(PARAMETERS(strcmp('R(Latm/Kmol)',PARAMETERS(:,1)),2));
Pm.V_l = str2double(PARAMETERS(strcmp('Vl',PARAMETERS(:,1)),2));
Pm.V_g = str2double(PARAMETERS(strcmp('Vg',PARAMETERS(:,1)),2));
Pm.dataset = str2double(PARAMETERS(strcmp('dataset',PARAMETERS(:,1)),2));
Pm.tstart = str2double(PARAMETERS(strcmp('tstart',PARAMETERS(:,1)),2));
Pm.tend = str2double(PARAMETERS(strcmp('tend',PARAMETERS(:,1)),2));
Pm.tstep = str2double(PARAMETERS(strcmp('tstep',PARAMETERS(:,1)),2));
Pm.ncompl = length(Pm.S);
Pm.ncompg = length(Pm.kla);
Pm.Htoti = find(strcmp('H+.tot',Comp.master));

% Maximum rate constants (corrected for temperature)
Pm.C_Ri = str2double(MT_l(3:end,3));
Pm.Kmax = str2double(MT_l(3:end,2));
k1 = find(strcmp('Yes', MT_l(3:end,4)));
Pm.Kmax(k1,1) = Pm.Kmax(k1,1).*exp(64000/8.314*(1/303-1/Pm.T)); % Veeken & Hamelers 1998

%  Initial total concentrations liquid & gas phase 
k1 = strfind(Comp.master,'.logact'); k1 = find(cellfun(@isempty,k1)==0);
Mout_add = cellfun(@(x) [x(1:end-7) '.tot'], Comp.master(k1), 'UniformOutput', false);
k2 = strfind(Comp.master,'pH'); k2 = find(cellfun(@isempty,k2)==0);
H_add = cellfun(@(x) 'H+.tot', Comp.master(k2), 'UniformOutput', false);
Mout_add = [Mout_add H_add]; k1 = [k1 k2]; Comp.all = [Comp.all Mout_add]; Comp.alli = [Comp.alli 0.001*ones(1,length(Mout_add))];
% Cini = initialize_ORI(Comp,1);
OrchestraObj = OrchestraCompositeCl();
OrchestraObj = OrchestraObj.add(OrchestraCl(Comp));
Cini = OrchestraObj.Calculate();

Comp.master(k1) = Mout_add; Comp.masteri(k1) = Cini{1}(end+1-length(k1):end);
Comp.all = [Comp.master Comp.out]; Comp.alli = [Comp.masteri Comp.outi];

MT_g_ini = Cini{1}(Pm.Ci_tf).*Pm.H.*Pm.V_g./Pm.R_gas./Pm.T; % mol px = Hx*Cx nx = pxVg/R
N2i = Pm.p*Pm.V_g/Pm.R_gas/Pm.T-sum(MT_g_ini(2:end));
MT_g_ini = [N2i;MT_g_ini(2:end)];
MT_l_ini = Comp.masteri*Pm.V_l;
Pm.names_CT = [Comp.master Pm.gas_name' 'CO2_{cum}' 'CH4_{cum}'];
Pm.names_CD = Comp.out ;

% Cumulative total concentration outputs 
k1 = strfind(Comp.out,'.cum'); k1 = find(cellfun(@isempty,k1)==0);
Mcum_add = cellfun(@(x) x(1:end-4), Comp.out(k1), 'UniformOutput', false);
Pm.cumgasi = find(ismember(Pm.gas_name,Mcum_add));  

% Initial masses ODE solver
IM = [MT_l_ini MT_g_ini' Comp.outi(k1)];

% Initialize ORCHESTRA
% ORI = initialize_ORI(Comp,0);
OrchestraObj = OrchestraCompositeCl;
OrchestraObj = OrchestraObj.add(OrchestraCl(Comp));

% Load experimental data
if Pm.dataset > 0; Meas = load_data(Pm.dataset); end

end

function [MT_l, MT_g, OUTPUT, PARAMETERS, INHIBITION_TYPE] = read_Pmatrix(fn)
 
ii = 0;
fid = fopen(fn); tline = fgetl(fid); 
while tline ~= -1,
    ii = ii+1;
    line = textscan(tline,'%s','delimiter',';');
    Pmatrix{ii,1} = line{:}';
    tline = fgetl(fid);
end
fclose (fid);

% Reads MT_l, total concentrations liquid phase
for j = 1:length(Pmatrix(:))
    
    if isempty(Pmatrix{j}{1}) == 1
        break
    end
    
    for i = 1:length(Pmatrix{j})
        if isempty(Pmatrix{j}{i}) == 1
            break
        else
            MT_l{j,i} = Pmatrix{j}{i};
        end
    end
    N = j;
end

% Reads output concentrations for ORCHESTRA
ii = 0; 
for j = N+1:length(Pmatrix(:))
    
    if isempty(Pmatrix{j}{1}) == 1 && ii > 0
        break
    elseif isempty(Pmatrix{j}{1}) == 1 
        N = N+1;
    end
    
    if isempty(Pmatrix{j}{1}) == 0
        ii = 1;
        for i = 1:length(Pmatrix{j})
            if isempty(Pmatrix{j}{i}) == 1
                break
            else
                OUTPUT{j-N,i} = Pmatrix{j}{i};
            end
        end
    end
end

% Reads PARAMETERS, process parameters
ii = 0; W = length(OUTPUT(1,:));
for j = N+1:length(Pmatrix(:))
    
    for i = 1+W:length(Pmatrix{j})
        
        if isempty(Pmatrix{j}{i}) == 0
            ii = 1;
            PARAMETERS{j-N,i-W} = Pmatrix{j}{i};
        end
        
        if isempty(Pmatrix{j}{i}) == 1 && ii == 0
            W = W+1;
        end
        
        if isempty(Pmatrix{j}{i}) == 1 && ii == 1
            break
        end
    end
end

% Reads INHIBITION_TYPE, inhibition constants and indices
ii = 0; W = W + length(PARAMETERS(1,:));
for j = N+1:length(Pmatrix(:))
    
    for i = 1+W:length(Pmatrix{j})
        
        if isempty(Pmatrix{j}{i}) == 0
            ii = 1;
            INHIBITION_TYPE{j-N,i-W} = Pmatrix{j}{i};
        end
        
        if isempty(Pmatrix{j}{i}) == 1 && ii == 0
            W = W+1;
        end
        
        if isempty(Pmatrix{j}{i}) == 1 && ii == 1
            break
        end
    end
end

% Reads MT_g, total concentrations gas phase 
ii = 0; W = W + length(INHIBITION_TYPE(1,:));
for j = N+1:length(Pmatrix(:))
    
    for i = 1+W:length(Pmatrix{j})
        
        if isempty(Pmatrix{j}{i}) == 0
            ii = 1;
            MT_g{j-N,i-W} = Pmatrix{j}{i};
        end
        
        if isempty(Pmatrix{j}{i}) == 1 && ii == 0
            W = W+1;
        end
        
        if isempty(Pmatrix{j}{i}) == 1 && ii == 1
            break
        end
    end
end   
end

function Meas = load_data(dataset)

% Load datafiles
data1 = load('../Data/file1'); data1 = data1.data; [data1{strcmp(data1,'null')}] = deal(NaN);
data2 = load('../Data/file2'); data2 = data2.data; [data2{strcmp(data2,'null')}] = deal(NaN);
data3 = load('../Data/file3'); data3 = data3.data; [data3{strcmp(data3,'null')}] = deal(NaN);

% Name data
Meas.Biogas_m = cell2mat(data1(2:end,9))./1000;  % m^3
Meas.EC_m = cell2mat(data1(2:end,5)); % mS/cm
Meas.CH4p_m = cell2mat(data1(2:end,10))/100; % percentage
Meas.CO2p_m = cell2mat(data1(2:end,11))/100; % percentage
Meas.pH_m = cell2mat(data1(2:end,4));
Meas.NH4_m = cell2mat(data2(2:end,4))/1000/18; % mol/L
Meas.VFA_m = cell2mat(data3(2:end,4))/1000/60.05; % mol/L

% Name times corresponding to datasets
tm = data1(2:end,1);  % t --> Biogas_m, pH_m, CH4p_m, CO2p_m, EC_m
tm2 = data2(2:end,1); % t --> NH4_m
tm3 = data3(2:end,1); % t --> VFA_m

switch dataset
    case 1
        Meas.Biogas_m = Meas.Biogas_m(1:224);
        Meas.EC_m = Meas.EC_m(1:224);
        Meas.CH4p_m = Meas.CH4p_m(1:224);
        Meas.CO2p_m = Meas.CO2p_m(1:224);
        Meas.pH_m = Meas.pH_m(1:224);
        Meas.NH4_m = Meas.NH4_m(2:132);
        Meas.VFA_m = Meas.VFA_m(1:56);
        tm = tm(1:224); tm = datenum(tm); fd = tm(1); tm = tm-fd+1;
        tm2 = tm2(2:132); tm2 = datenum(tm2); tm2 = tm2-fd+1;
        tm3 = tm3(1:56); tm3 = datenum(tm3); tm3 = tm3-fd+1;
        
    case 2
        Meas.Biogas_m = Meas.Biogas_m(226:386);
        Meas.EC_m = Meas.EC_m(226:386);
        Meas.CH4p_m = Meas.CH4p_m(226:386);
        Meas.CO2p_m = Meas.CO2p_m(226:386);
        Meas.pH_m = Meas.pH_m(226:386);
        Meas.NH4_m = Meas.NH4_m(133:234);
        Meas.VFA_m = Meas.VFA_m(58:100);
        tm = tm(226:386); tm = datenum(tm); fd = tm(1); tm = tm-fd+1;
        tm2 = tm2(133:234); tm2 = datenum(tm2); tm2 = tm2-fd+1;
        tm3 = tm3(58:100); tm3 = datenum(tm3); tm3 = tm3-fd+1;
        
    case 3
        Meas.Biogas_m = Meas.Biogas_m(387:464);
        Meas.EC_m = Meas.EC_m(387:464);
        Meas.CH4p_m = Meas.CH4p_m(387:464);
        Meas.CO2p_m = Meas.CO2p_m(387:464);
        Meas.pH_m = Meas.pH_m(387:464);
        Meas.NH4_m = Meas.NH4_m(235:277);
        Meas.VFA_m = Meas.VFA_m(102:131);
        tm = tm(387:464); tm = datenum(tm); fd = tm(1); tm = tm-fd+1;
        tm2 = tm2(235:277); tm2 = datenum(tm2); tm2 = tm2-fd+1;
        tm3 = tm3(102:131); tm3 = datenum(tm3); tm3 = tm3-fd+1;
        
    case 4
        Meas.Biogas_m = Meas.Biogas_m(465:542);
        Meas.EC_m = Meas.EC_m(465:542);
        Meas.CH4p_m = Meas.CH4p_m(465:542);
        Meas.CO2p_m = Meas.CO2p_m(465:542);
        Meas.pH_m = Meas.pH_m(465:542);
        Meas.NH4_m = Meas.NH4_m(278:320);
        Meas.VFA_m = Meas.VFA_m(132:161);
        tm = tm(387:464); tm = datenum(tm); fd = tm(1); tm = tm-fd+1;
        tm2 = tm2(278:320); tm2 = datenum(tm2); tm2 = tm2-fd+1;
        tm3 = tm3(132:161); tm3 = datenum(tm3); tm3 = tm3-fd+1;
        
    case 5
        Meas.Biogas_m = Meas.Biogas_m(543:620);
        Meas.EC_m = Meas.EC_m(543:620);
        Meas.CH4p_m = Meas.CH4p_m(543:620);
        Meas.CO2p_m = Meas.CO2p_m(543:620);
        Meas.pH_m = Meas.pH_m(543:620);
        Meas.NH4_m = Meas.NH4_m(321:362);   % skipped last data point (three months later than all other)
        Meas.VFA_m = Meas.VFA_m(162:191);
        tm = tm(543:620); tm = datenum(tm); fd = tm(1); tm = tm-fd+1;
        tm2 = tm2(321:362); tm2 = datenum(tm2); tm2 = tm2-fd+1;
        tm3 = tm3(162:191); tm3 = datenum(tm3); tm3 = tm3-fd+1;
        
    case 6
        Meas.Biogas_m = Meas.Biogas_m(621:698);
        Meas.EC_m = Meas.EC_m(621:698);
        Meas.CH4p_m = Meas.CH4p_m(621:698);
        Meas.CO2p_m = Meas.CO2p_m(621:698);
        Meas.pH_m = Meas.pH_m(621:698);
        Meas.NH4_m = Meas.NH4_m(364:405); % skipped last data point (three months later than all other)
        Meas.VFA_m = Meas.VFA_m(192:221);
        tm = tm(621:698); tm = datenum(tm); fd = tm(1); tm = tm-fd+1;
        tm2 = tm2(364:405); tm2 = datenum(tm2); tm2 = tm2-fd+1;
        tm3 = tm3(192:221); tm3 = datenum(tm3); tm3 = tm3-fd+1;
        
    case 7
        Meas.Biogas_m = Meas.Biogas_m(699:776);
        Meas.EC_m = Meas.EC_m(699:776);
        Meas.CH4p_m = Meas.CH4p_m(699:776);
        Meas.CO2p_m = Meas.CO2p_m(699:776);
        Meas.pH_m = Meas.pH_m(699:776);
        Meas.NH4_m = Meas.NH4_m(407:449);
        Meas.VFA_m = Meas.VFA_m(222:251);
        tm = tm(699:776); tm = datenum(tm); fd = tm(1); tm = tm-fd+1;
        tm2 = tm2(407:449); tm2 = datenum(tm2); tm2 = tm2-fd+1;
        tm3 = tm3(222:251); tm3 = datenum(tm3); tm3 = tm3-fd+1;
        
    case 8
        Meas.Biogas_m = Meas.Biogas_m(777:854);
        Meas.EC_m = Meas.EC_m(777:854);
        Meas.CH4p_m = Meas.CH4p_m(777:854);
        Meas.CO2p_m = Meas.CO2p_m(777:854);
        Meas.pH_m = Meas.pH_m(777:854);
        Meas.NH4_m = Meas.NH4_m(450:492);
        Meas.VFA_m = Meas.VFA_m(252:281);
        tm = tm(777:854); tm = datenum(tm); fd = tm(1); tm = tm-fd+1;
        tm2 = tm2(450:492); tm2 = datenum(tm2); tm2 = tm2-fd+1;
        tm3 = tm3(252:281); tm3 = datenum(tm3); tm3 = tm3-fd+1;
        
    case 9
        Meas.Biogas_m = Meas.Biogas_m(855:932);
        Meas.EC_m = Meas.EC_m(855:932);
        Meas.CH4p_m = Meas.CH4p_m(855:932);
        Meas.CO2p_m = Meas.CO2p_m(855:932);
        Meas.pH_m = Meas.pH_m(855:932);
        Meas.NH4_m = Meas.NH4_m(493:534); % skipped last data point (three months later than all other)
        Meas.VFA_m = Meas.VFA_m(282:311);
        tm = tm(855:932); tm = datenum(tm); fd = tm(1); tm = tm-fd+1;
        tm2 = tm2(493:534); tm2 = datenum(tm2); tm2 = tm2-fd+1;
        tm3 = tm3(282:311); tm3 = datenum(tm3); tm3 = tm3-fd+1;
end

% Remove all the NaN in the data & make datasets unique
Meas.Biogas_tm = tm;   Meas.Biogas_tm(isnan(Meas.Biogas_m)) = [];   Meas.Biogas_m(isnan(Meas.Biogas_m)) = [];
[Meas.Biogas_tm, id] = unique(Meas.Biogas_tm);  Meas.Biogas_m = Meas.Biogas_m(id);
Meas.pH_tm = tm;    Meas.pH_tm(isnan(Meas.pH_m)) = [];     Meas.pH_m(isnan(Meas.pH_m)) = [];
[Meas.pH_tm, id] = unique(Meas.pH_tm);  Meas.pH_m = Meas.pH_m(id);
Meas.CH4p_tm = tm;  Meas.CH4p_tm(isnan(Meas.CH4p_m)) = []; Meas.CH4p_m(isnan(Meas.CH4p_m)) = [];
[Meas.CH4p_tm, id] = unique(Meas.CH4p_tm);  Meas.CH4p_m = Meas.CH4p_m(id);
Meas.CO2p_tm = tm;  Meas.CO2p_tm(isnan(Meas.CO2p_m)) = []; Meas.CO2p_m(isnan(Meas.CO2p_m)) = [];        
[Meas.CO2p_tm, id] = unique(Meas.CO2p_tm);  Meas.CO2p_m = Meas.CO2p_m(id);
Meas.NH4_tm = tm2;  Meas.NH4_tm(isnan(Meas.NH4_m)) = [];   Meas.NH4_m(isnan(Meas.NH4_m)) = [];
[Meas.NH4_tm, id] = unique(Meas.NH4_tm);  Meas.NH4_m = Meas.NH4_m(id);
Meas.VFA_tm = tm3;  Meas.VFA_tm(isnan(Meas.VFA_m)) = [];   Meas.VFA_m(isnan(Meas.VFA_m)) = [];
[Meas.VFA_tm, id] = unique(Meas.VFA_tm);  Meas.VFA_m = Meas.VFA_m(id);
end