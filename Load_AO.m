%% Load AlphaOmega data
%
% Subfunction for Anodic_Cathodic_main
%
% Loads raw data from .mat files saved with the AlphaOmega SnR system
%
function data = Load_AO(filenames,noPar,varargin)
% Function to load data recorded with the AlphaOmega AlphaLabSnR-System and
% converted to the Matlab ".mat" file format, using the official AlphaOmega
% converter.
%
% Usage:
%       data = Load_AO(filenames,noPar,'option','value')
%
% filenames = String or cell array of strings containing the filename(s)
% to be loaded
%
% noPar = Number of different parameter levels to be loaded (e.g.
% intensities of stimulation)
%
% Additionally the following Option/Value pairs may be provided
% (case-sensitive):
%
% probe = Number of shank or probe to be loaded.
%         probe = 1 -> Channels  1-16
%         probe = 2 -> Channels 17-32
%         probe = 3 -> Channels 33-48
%         probe = 4 -> Channels 49-64
%         probe = NaN -> Load Current Monitor channel
%
% sweepBefore & sweepAfter = Number of samples to load before and after
% stimulation onset (TTL pulse or StimMarker, Default: 50 & 150 [in ms])
%
% sweeps = Number of stimulus repetitions per parameter level
% (Default: 32)
%
% tsAdjust = Adjust the timestamps vector by this vector
%            (e.g. tsAdjust = [1:50] only select timestamps 1 to 50)
%
% verbose = [0 = "no" | 1 = "yes"] Verbose messages about
% loading progress (Default: 0)
%

%Check inputarguments
options = struct(...
    'sweepBefore',50,...
    'sweepAfter',150,...
    'sweeps',32,...
    'probe',1,...
    'channels',[],...
    'tsAdjust',[],...
    'loadTS',0,...
    'verbose',0);

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

% Count number of files provided
if ~iscell(filenames)
    noOfFiles = 1;
    filetmp{1} = filenames;
    filenames = filetmp;
elseif iscell(filenames) && numel(filenames) == 1
    noOfFiles = 1;
else
    noOfFiles = numel(filenames);
end

% Determine channels to load
if isempty(options.channels)
    if isnan(options.probe)
        options.channels = 1;
    else
        options.channels = (16*options.probe)-15:(16*options.probe);
    end
end
noOfChannels = length(options.channels);

% First electrode
firstElectrode = sprintf('%02d',options.channels(1));

if options.verbose
    fprintf('Loading file 01 of %02.0f\n',noOfFiles)
end

% Number of sweeps supposed to load:
supp_sweeps = options.sweeps*noPar;

% Initialise variables
for ch = 1:noOfChannels
%     v = sprintf('CSPK_0%02d',options.channels(ch));
    fileData.(['ch' num2str(ch)]) = [];
end

% Loop through all files and put data together to one long recording
for fileIDX = 1:noOfFiles
    
    % Load data from file
    load(filenames{fileIDX});
    
    if exist('CStimMarker_001','var')== 1
        eval(['tmstmps_temp{fileIDX} = (CStimMarker_001(1,:)/(CStimMarker_001_KHz*1000))-CSPK_0' firstElectrode '_TimeBegin;']);
    elseif exist('CTTL_001_Up','var')== 1
        eval(['timeoffset{fileIDX} = (CTTL_001_TimeBegin)-(CSPK_0' firstElectrode '_TimeBegin);']);
        tmstmps_temp{fileIDX} = (CTTL_001_Up(1,:)/(CTTL_001_KHz*1000))+timeoffset{fileIDX}; %#ok
    else
        if options.verbose
            fprintf('Skipping file %s (File without stimulus!)\n',filenames{fileIDX})
        end
        noOfFiles = noOfFiles-1;
        continue
    end
    
       % Current monitor channel handling
    if isnan(options.probe)
        if exist('CSPK_049','var')
            options.channels = 49;
            firstElectrode = '49';
        elseif exist('CSPK_033','var')
            options.channels = 33;
            firstElectrode = '33';
        else
            options.channels = 17;
            firstElectrode = '17';
        end
    end
    
    eval(['tm{fileIDX} = round(tmstmps_temp{fileIDX}*(CSPK_0' firstElectrode '_KHz*1000));']);
    eval(['offset_temp{fileIDX} = length(CSPK_0' firstElectrode ');']);
    
    eval(['sweepBefore = (options.sweepBefore*(CSPK_0' firstElectrode '_KHz));']);
    eval(['sweepAfter = (options.sweepAfter*(CSPK_0' firstElectrode '_KHz));']);
    noSamples = sweepBefore+sweepAfter+1;
    
    for ch = 1:noOfChannels
        if options.verbose
            fprintf('File: %02d | %02d of %02d.\n',fileIDX,ch,noOfChannels);
        end
        v =  sprintf('CSPK_0%02d',options.channels(ch));
        fileData.(['ch' num2str(ch)]) = [fileData.(['ch' num2str(ch)]) int16(eval(v))];
    end
    clear('-regexp','^CSPK|^CStim|^CTTL');
    
    if options.verbose && fileIDX < noOfFiles
        fprintf('Loading file %2d of %2d\n',fileIDX+1,noOfFiles)
    end
end

% Concatenate timeStamps
timeStamps = tm{1}; %#ok
offset=cumsum(cell2mat(offset_temp));
for fileIDX = 2:noOfFiles
    timeStamps = [timeStamps tm{fileIDX}+offset(fileIDX-1)];  %#ok
end
clear tm tmstmps_temp timeoffset offset v ch

% Handling of wrong number of timestamps
if supp_sweeps ~= length(timeStamps)
    fprintf('ATTENTION! Wrong number of timestamps provided! Resulting data might be corrupted!\n')
    if options.verbose
        fprintf('ATTENTION! %d timestamps assumed (= %d parameters x %d repetitions), but %d timestamps found!\n',supp_sweeps,noPar,options.sweeps,length(timeStamps));
    end
    if ~isempty(options.tsAdjust)
        fprintf('ATTENTION! Timestamps adjustment vector found. Correcting timestamps...\n')
        eval(['timeStamps = timeStamps(' options.tsAdjust ');']);
    end
end

if ~options.loadTS
    if options.verbose
        fprintf('Attempting to create %d x %d x %d x %d matrix\n',noPar,noOfChannels,options.sweeps,noSamples)
    end
    data = zeros(noPar,noOfChannels,options.sweeps,noSamples,'int16');
    % Reshape to different strengths
    for str = 1:noPar
        for ch = 1:noOfChannels
            %         v =  sprintf('CSPK_0%02d',options.channels(ch)); %#ok
            for sw = 1:options.sweeps
                data(str,ch,sw,:) = eval(['fileData.([''ch'' num2str(ch)])(timeStamps((str*' num2str(options.sweeps) ')-' num2str(options.sweeps) '+sw)-sweepBefore:timeStamps((str*' num2str(options.sweeps) ')-' num2str(options.sweeps) '+sw)+sweepAfter)']);
            end
        end
    end
    clear fileData
else % load all timestamps and ignore par&sweeps
    if options.verbose
        fprintf('Attempting to create 1 x %d x %d x %d matrix\n',1,noOfChannels,length(timestamps),noSamples)
    end
    data = zeros(1,noOfChannels,length(timestamps),noSamples,'int16');
    for ch = 1:noOfChannels
        for sw = 1:length(timestamps)
            data(1,ch,sw,:) = fileData.(['ch' num2str(ch)])(timeStamps(sw)-sweepBefore:timeStamps(sw)+sweepAfter);
        end
    end
    clear fileData
end
    
    
% This is to make really sure, that the first and last sweep is removed from
% electrically stimulated data, which has usually 32 repetitions...
if size(data,3) == 32
    data = data(:,:,2:31,:);
end

if options.verbose
    fprintf('Data succesfully loaded!\n')
end
end