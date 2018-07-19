%% SPIKES
%
% Subfunction for Anodic_Cathodic_main/Load_MUA_Electric_VarDepth
%
% Detects spiking events as threshold crossings
%
% Written with significant input from Dr. Peter Hubka, Hannover Medical School
%
function [spikematrixlog] = spikes(data,threshold)

% Use default threshold value if none is given
if nargin < 2
    threshold = 4;
end

% Determine the size of the available data
[m,n,p] = size(data);

% Threshold criteria
% after Quiroga, 2004
thr = threshold*(median(abs(data/0.6745),3));

% Spike detection

% Initialize variable
spikematrix = zeros(m,n,p);

for c = 1:m
    for s = 1:n
        spkind = find(data(c,s,:)<-thr(c,s)); % find threshold crossings
        spkindhelper = [0; find(diff(spkind)>1); length(spkind)]; % only look at non-continuous 'spikes'
        
        if ~isempty(spkind) % make sure there are spikes
            
            % find maximal value for each spike
            spkindh = zeros(length(spkindhelper)-1,1);
            for i = 1:length(spkindhelper)-1
                spkindh(i)= spkindhelper(i)+find(data(c,s,spkind(spkindhelper(i)+1):spkind(spkindhelper(i+1)))==-max(abs(data(c,s,spkind(spkindhelper(i)+1):spkind(spkindhelper(i+1))))));
            end
            
            % write maximal value into spikematrix
            for i = 1:length(spkindh)
                spikematrix(c,s,spkind(spkindh(i))) = data(c,s,spkind(spkindh(i)));
            end
            
        end %if
        
        clear spkind spkindhelper spkindh
    end
end

% logical spike matrix: 1 at timepoint with spike, 0 elsewhere
spikematrixlog = spikematrix ~= 0;
