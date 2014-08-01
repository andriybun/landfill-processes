% Modify text file containing atmospheric boundary condition
function ModifyAtmosphIn(Extra)


% Open ATMOSPH.IN
fid = fopen([Extra.subdir '\H1D\ATMOSPH.IN'],'r+');

% Go to atmospheric information section
flag = [];
while isempty(flag)
	str = fgetl(fid);
	flag = findstr(str,'tAtm');
end

% Insert new lower boundary condition
fseek(fid,0,'cof');
fprintf(fid,'%11.0f %11.2f %11.4f %11.0f %11.0f %11.0f %11.2f %11.0f\n',Extra.boundcon');
fseek(fid,0,'eof');

% Close file
fclose(fid);

