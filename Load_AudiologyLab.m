%% Load AudiologyLab data
%
% Subfunction for Anodic_Cathodic_main
%
% Loads data files saved by AudiologyLab (Otoconsult)
%
function data = Load_AudiologyLab(file,chNumber,header)

% Number of channels to load
if nargin < 2
    chNumber = 16;
end

% Load only header information and no single sweep data
if nargin < 3
    header = 0;
end

% Open experimental datafile and extract relevant info

fid = fopen(file);

data.version = textscan(fid, '%*s%n',1);
data.version = data.version{:};
data.filename = textscan(fid, '%*s%s',1);
data.filename = data.filename{:};
textscan(fid,'%*s',25);
data.invert = textscan(fid, '%n',1);
data.invert = data.invert{:};
data.NoOfSweep = textscan(fid, '%*s%n',1);
data.NoOfSweep = data.NoOfSweep{:};
textscan(fid,'%*s',7);
data.SampleRate = textscan(fid, '%n',1);
data.SampleRate = data.SampleRate{:};       % in [kHz]!
textscan(fid,'%*s',12);
data.Duration = textscan(fid, '%n',1);
data.Duration = data.Duration{:};

fclose(fid);

% If data is collected with inverting stimuli the total number of stimulus
% repetitions is twice the number of repetitions set in the software
if data.invert == 1
    data.NoOfSweep = data.NoOfSweep * 2;
end

% Open parameter file and extract relevant information
parafile = strcat(file,'para');
fid = fopen(parafile);

% Things to ignore
textscan(fid, '%*n',36);
textscan(fid, '%*s',2);
textscan(fid, '%*n',36);
textscan(fid, '%*s',2);
textscan(fid, '%*n',76);

data.ChNumber = textscan(fid, '%n',1);
data.ChNumber = data.ChNumber{:};
fclose(fid);

% Open stimulation parameters file to sort single sweeps
% accordingly
stimparafile = strcat(file,'stimpara');

data.stimpara = dlmread(stimparafile,';');

% Sort stimulation parameters according to:
% Frequencies
[FreqLvl, ~, ~] = unique(data.stimpara(:,3));
% Attenuations
[AttLvl, ~, ~] = unique(data.stimpara(:,4));

% Open single sweep data file to extract raw data
if header == 0
sweepsfile = strcat(file,'sweeps');

fid = fopen(sweepsfile);

for i = 1:length(FreqLvl)
    for j = 1:length(AttLvl)
        for k = 1:data.NoOfSweep
            for l = 1:chNumber
                data.data(:,l,k,j,i) = fread(fid,(data.Duration/1000)*(data.SampleRate*1000),'double=>single');
            end
        end
    end
end
end

end