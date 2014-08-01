% Modify text file containing soil hydraulic parameters
function ModifySelectorIn(x,Extra)


% Open SELECTOR.IN
fid = fopen([Extra.subdir '\H1D\SELECTOR.IN'],'r+');

% Go to water flow section
flag = [];
while isempty(flag)
	str = fgetl(fid);
	flag = findstr(str,'thr');
end

% Insert water flow parameters
fseek(fid,0,'cof');
fprintf(fid,' %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e\r\n',x(1:6));
fseek(fid,0,'eof');

% Close file
fclose(fid);


