%% Calculate current-source density
%
% Subfunction for Anodic_Cathodic_main
%
function exps_out = findexperiments(varargin)
%
% Searches in the folder *searchFolder* for the presence of files matching 
% *expression*. Returns a cell array of strings containing experiment names
% 
% *searchFolder* defaults to: 'K:\ICMS-GP\Data\GP*'
% *expression* defaults to: '*'
%


%Check inputarguments
options = struct(...
    'expression','*',...
    'searchFolder','K:\ICMS-GP\Data\GP*');

optionNames = fieldnames(options);
nOptionalArgs = length(varargin);
if round(nOptionalArgs/2)~=nOptionalArgs/2
    error('Please input Name/Value pairs!');
end

for pair = reshape(varargin,2,[])
    inpName = pair{1}; 
    if any(strcmp(inpName,optionNames))
        options.(inpName) = pair{2};
    else
        error('%s is not a valid parameter name',inpName)
    end
end

exps_out = {};
expFolder = dir(options.searchFolder);

for idx = 1:length(expFolder)
    tempFolder = dir(['K:\ICMS-GP\Data\' expFolder(idx).name '\AlphaOmega\MATFiles\' options.expression]);
    if ~isempty(tempFolder)
        exps_out = [exps_out expFolder(idx).name]; %#ok don't know beforehand how many will be found
    end
end
