% Read text file containing the time series of simulated soil water contents
function sims = ReadObsNodeOut(Extra)

% Open OBS_NODE.OUT
fid = fopen([Extra.subdir '\H1D\OBS_NODE.OUT']);

% Go to data section
flag = [];
while isempty(flag)
	str = fgetl(fid);
	flag = findstr(str,'time');
end

% Read simulated soil water contents
cols = 4;
rows = Inf;
data = fscanf(fid,'%f',[cols rows]);

% Store simulated soil water contents in structure
data = data';
sims.hoy = data(:,1);
sims.water = data(:,3);

% close file
fclose(fid);

