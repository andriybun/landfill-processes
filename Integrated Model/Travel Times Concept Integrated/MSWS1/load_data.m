function Meas = load_data(Pm)

R = Pm(3); T = Pm(1); p = Pm(2); 

% load experimental data in XX2
a = fopen('Data_MSWS1.csv');
N = 14; M = 225; XX2 = [];
for i=1:M
    X = textscan(a,'%s',N);
    XX2 = [XX2;X{:}'];
end
fclose('all');

% Name data
Meas.Bio_m = str2double(XX2(2:225,2))/(R*T/p)/1000;
Meas.pH_m = str2double(XX2(2:225,5));
Meas.VFA_m = str2double(XX2(2:58,7))/1000/60.05;
Meas.NH4_m = str2double(XX2(2:132,9))/1000/18;  
Meas.EC_m = str2double(XX2(2:225,12)); % mS/cm
Meas.CH4p_m = str2double(XX2(2:225,13))/100; % percentage
Meas.CO2p_m = str2double(XX2(2:225,14))/100; % percentage

tm = str2double(XX2(2:225,1));  % t --> Biom, pHm, CH4p_m, CO2p_m, EC
tm2 = str2double(XX2(2:58,6)); % t --> VFAm
tm3 = str2double(XX2(2:132,8)); % t --> NH4m

% Remove all the NaN in the data & strange values in EC
Meas.Bio_tm = tm;   Meas.Bio_tm(isnan(Meas.Bio_m)) = [];   Meas.Bio_m(isnan(Meas.Bio_m)) = []; 
Meas.pH_tm = tm;    Meas.pH_tm(isnan(Meas.pH_m)) = [];     Meas.pH_m(isnan(Meas.pH_m)) = [];
Meas.VFA_tm = tm2;  Meas.VFA_tm(isnan(Meas.VFA_m)) = [];   Meas.VFA_m(isnan(Meas.VFA_m)) = [];
Meas.NH4_tm = tm3;  Meas.NH4_tm(isnan(Meas.NH4_m)) = [];   Meas.NH4_m(isnan(Meas.NH4_m)) = [];
Meas.EC_tm = tm;    Meas.EC_tm(isnan(Meas.EC_m)) = [];     Meas.EC_m(isnan(Meas.EC_m)) = [];
Meas.CH4p_tm = tm;  Meas.CH4p_tm(isnan(Meas.CH4p_m)) = []; Meas.CH4p_m(isnan(Meas.CH4p_m)) = [];
Meas.CO2p_tm = tm;  Meas.CO2p_tm(isnan(Meas.CO2p_m)) = []; Meas.CO2p_m(isnan(Meas.CO2p_m)) = [];
end