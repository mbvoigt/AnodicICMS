%% mbv - Custom logging function
function [ret] = mbv_log(msg,varargin)

%Check inputarguments
options = struct('file',[],...
    'type','   ',...
    'time',1);

optionNames = fieldnames(options);
nOptionalArgs = length(varargin);
if round(nOptionalArgs/2)~=nOptionalArgs/2
    error('Please input Name/Value pairs!');
end

for pair = reshape(varargin,2,[])
    inpName = pair{1}; 
    if any(strcmpi(inpName,optionNames))
        options.(inpName) = pair{2};
    else
        error('%s is not a valid parameter name',inpName)
    end
end

% Convert type code to message
if options.type == 0
    options.type = 'ERR';
elseif options.type == 1
    options.type = '   ';
elseif options.type == 2
    options.type = 'WRN';
end

% Log message
if isempty(options.file) % display in command window
    if options.time
        ret = ['[' datestr(now,'yyyy-mm-dd-HH:MM:SS') '] ' options.type ' | ' msg];
        fprintf('%s\n',ret) 
    else
        ret = msg;
        fprintf('%s\n',ret)
    end
else % write log file
    fID = fopen(options.file,'at');
    if options.time
        ret = ['[' datestr(now,'yyyy-mm-dd-HH:MM:SS') '] ' options.type ' | ' msg];
        fprintf(fID,'%s\n',ret);
    else
        ret = msg;
        fprintf(fID,'%s\n',ret);
    end
    fclose(fID);
end
    

