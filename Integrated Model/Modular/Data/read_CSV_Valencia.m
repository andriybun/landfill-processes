function [data] = read_CSV_Valencia(fn)

addpath(genpath('.'));
fid = fopen(fn);

%read header line
% tline = fgetl(fid); % skip first line
tline = fgetl(fid); % header line

%parse to header cell matrix
headers_rd = textscan(tline,'%s','delimiter',';');
headers = headers_rd{1};
% headers = ['date'; 'time'; headers(2:end)]; % split date and time
ntype = length(headers);

%read next line from file until eof
ii = 0;
tline = fgetl(fid);
while tline ~= -1,
    ii = ii+1;
    datastr(ii) = textscan(tline,'%s','delimiter',';');
    tline = fgetl(fid);
end
fclose (fid);
ndata = ii;

data = cell(ndata,ntype);
[data{:}] = deal('null');
for ii = 1:ndata,
%     datetime = datastr{ii}{1};  date = datetime(1:10); date = [date(7:end) date(3:5) date(6) date(1:2)]; time = datetime(12:end);
%     data{ii,1} = date; data{ii,2} = time;
    for jj = [1 2 length(datastr{ii})]    
        input = datastr{ii}{jj};
        if isempty(input) ~= 1
           data{ii,jj} = input;
        end
    end
    for jj = 3:length(datastr{ii})-1
        input = str2double(datastr{ii}{jj});
        if isnan(input) ~= 1
            data{ii,jj} = input;
        end
    end
end

data = [headers';data];

save(fn(1:end-4),'data'); 
end
