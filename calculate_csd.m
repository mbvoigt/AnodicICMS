%% Calculate current-source density
%
% Subfunction for Anodic_Cathodic_main
%
% Calculates current-source densities from local field potentials
%
function [csd,varargout] = calculate_csd(data,stepSize,elStimRemove,stimCh,probe)

% Distance between contacts
if nargin < 2
    stepSize = 0.15;
end

% If set to one: remove stimulation electrode
if nargin < 3
    elStimRemove = 0;
end

% Channel number to remove
if nargin < 4
    stimCh = [];
end

% Use 16 channel or 32 channel electrode map
if nargin < 5
    probe = 16;
end

if elStimRemove == 1 && isempty(stimCh)
    error('Stimulation channel has to be specified in order to be removed!')
end

if probe == 16
    channelidx = [9,8,10,7,13,4,12,5,15,2,16,1,14,3,11,6];
elseif probe == 32
    channelidx = [9,8,10,7,11,6,12,5,13,4,14,3,15,2,16,1];
else
   channelidx = [9,8,10,7,13,4,12,5,15,2,16,1,14,3,11,6]; 
end

% Reorder channels according to electrode channel map
data = data(channelidx,:);

% Add Vaknin electrodes
datatmp = data;
data(1,:) = datatmp(1,:);
for t = 2:17
    data(t,:) = datatmp(t-1,:);
end
data(18,:) = datatmp(16,:);


[m,n] = size(data);
csd = zeros(m-2,n);

if elStimRemove
    %calculate stepSize & remove stimchannels
    stSize = repmat(stepSize,m,1);
    k=1;
    for ch = 1:m
        if ~ismember(ch,stimCh)
            data_temp(k,:) = data(ch,:);
            k=k+1;
        else
            if ch > 1
                stSize(ch-1)=stSize(ch-1)+stepSize;
            end
        end
    end
    %calculate csd_temp
    csd_temp = zeros(m-length(stimCh)-2,n);
    for j = 2:size(data_temp,1)-1
        csd_temp(j-1,:) = (data_temp(j+1,:)-(2*data_temp(j,:))+data_temp(j-1,:))/((stSize(j))^2);
    end

    %insert stimchannel 0s
    if ismember(1,stimCh) && numel(stimCh)==1
        for ch = 2:m-2
            csd(ch,:) = csd_temp(ch-1,:);
        end
    elseif ismember(16,stimCh) && numel(stimCh)==1
        for ch = 1:size(csd_temp,1)
            csd(ch,:) = csd_temp(ch,:);
        end
    else
        k=1;
        for ch = 1:size(csd,1)
            if k > size(csd_temp,1)
                continue
            else
                if ismember(ch,stimCh)
                    continue
                elseif ~ismember(ch,stimCh)
                    csd(ch,:) = csd_temp(k,:);
                    k=k+1;
                else
                    continue
                end
            end
        end
    end
else
    for j = 2:m-1
        csd(j-1,:) = (data(j+1,:)-(2*data(j,:))+data(j-1,:))/(stepSize^2);
    end

end

%%%% This is to invert the CSD!
% Afterwards sinks should be positive and sources should be negative!
csd = -csd;
%%%%

if nargout > 1 % include AVREC and residuals
    varargout(1) = {sum(abs(csd),1)/14};
    varargout(2) = {(sum(csd,1))./(sum(abs(csd),1))};
end