%% BLANKSTIMULUS
%
% Subfunction for Anodic_Cathodic_main/Load_CSD_Electric_IO_Electrode_[...]
%
% Removes stimulation artefact by linear interpolation in a given range
%
function blankedsignal = blankstimulus(data,SR,blankstart,blankstop)
% Function to blank the stimulus artefact in neural recordings using a linear interpolation method.
%
% USAGE:
%
%   Blankedsignal = blankstimulus(Data,SamplingRate, Start, Stop)
%
% "Data" should be the standard data format of the ICMS experiments,
% meaning a matrix in the form of data(channels,sweeps,"rawdata" == voltages)
% "SamplingRate" in Hz 
% "Start" in ms
% "Stop" in ms

blankstart = blankstart/1000;
blankstop = blankstop/1000;

blfrom = ceil(blankstart*SR); 
blto = ceil(blankstop*SR);
bldur = blto - blfrom+1;

blankedsignal = data;

try
[channel, sweep, ~] = size(data);

for c = 1:channel
    for s = 1:sweep
        blankedsignal(c,s,blfrom:blto) = ...
            linspace(double(data(c,s,blfrom)),double(data(c,s,blto)),...
            bldur);
    end
end

catch exception
    disp('Data format does not match the expected ICMS format.')
    disp(exception)
end