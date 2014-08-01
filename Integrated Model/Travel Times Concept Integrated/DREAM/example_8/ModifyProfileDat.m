% Modify text file containing initial condition
function ModifyProfileDat(Extra)


% Open PROFILE.DAT
fid = fopen([Extra.subdir '\H1D\PROFILE.DAT'],'r+');

% Go to profile section
flag = [];
while isempty(flag)
	str = fgetl(fid);
	flag = findstr(str,'Mat');
end

% Insert new initial condition
fseek(fid,0,'cof');
fprintf(fid,'%5.0f %15.6e %15.6e %3.0f %3.0f %15.6e %15.6e %15.6e %15.6e\n',Extra.initial');
fseek(fid,0,'eof');

% Close file
fclose(fid);