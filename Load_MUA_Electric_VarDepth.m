%% Load MUA electric - Varying stimulation depth
%
% Subfunction for Anodic_Cathodic_main
%
% Loads multi-unit activity data for electrically stimulated trials with
% varying stimulation depth
%
function [MUAA,MUAC,MUASpikesA,MUASpikesC] = Load_MUA_Electric_VarDepth(cfg)

% Define experiments with complete anodic- and cathodic-leading dataset
cfg.Anoexps = cfg.Anoexps([1,3,5,6,7,8,10,11]);

% Define Neuronexus channel mapping for experiment which used the 1x16
% probe instead of the 2x16 probes
tmpchannelidx = [9,8,10,7,13,4,12,5,15,2,16,1,14,3,11,6];

% Preallocation
MUAA = zeros(length(cfg.Anoexps),16,16,5501);
MUAC = zeros(length(cfg.Anoexps),16,16,5501);
MUASpikesA = zeros(length(cfg.Anoexps),16,16,30,5501);
MUASpikesC = zeros(length(cfg.Anoexps),16,16,30,5501);

% Repeat for every experiment
for exps = 1:length(cfg.Anoexps)
    
    % Define folder where data can be found for this specific experiment
    fprintf('%2d | 11 - %s | anodic',exps,cfg.Anoexps{exps})
    pth = ['K:\ICMS-GP\Data\' cfg.Anoexps{exps} '\AlphaOmega\MATFiles\'];
    position = 1; % if more than one position would have been measured this could be implemented here
    
    %----------------------------
    % Anodic-leading stimulation
    %----------------------------
    
    % Determine filenames according to experiment to load and load raw data
    data = zeros(16,16,30,5501);
    for el = 1:16
        if strcmp(cfg.Anoexps{exps},'GP09F16CL') || strcmp(cfg.Anoexps{exps},'GP13I16CL')
            file = sprintf('%sctx_pos_%d_layering_100%02d_anod0001.mat',pth,position,cfg.channelidx(el)-1);
        elseif strcmp(cfg.Anoexps{exps},'GP11K15CL')
            file = sprintf('%sctx_pos%d_layering_anodic_100%02d0001.mat',pth,position,tmpchannelidx(el)-1);
        else
            file = sprintf('%sctx_pos%d_layering_100%02d_anodi0001.mat',pth,position,cfg.channelidx(el)-1);
        end
        data(el,:,:,:) = Load_AO(file,1,'sweepBefore',50,'sweepAfter',200,'sweeps',32);
    end
    
    % Stimulation artefact blanking
    dataBlanked = zeros(16,16,30,5501);
    for el = 1:16
        dataBlanked(el,:,:,:) = ...
            blankstimulus(squeeze(data(el,:,:,:)),22000,50,53);
    end
    clear data
    
    % Filter for MUA
    for el = 1:size(dataBlanked,1)
        for ch = 1:size(dataBlanked,2)
            for sw = 1:size(dataBlanked,3)
                data(el,ch,sw,:) = filtfilt(cfg.filter.bS,cfg.filter.aS,double(squeeze(dataBlanked(el,ch,sw,:))));
            end
        end
    end
    clear el ch sw dataBlanked
    
    % Determine spikes
    for el = 1:16
        [~,MUASpikesA(exps,el,:,:,:)] = spikes(squeeze(data(el,:,:,:)),4);
    end
    MUAA(exps,:,:,:,:) = data;
    clear data

    %------------------------------
    % Cathodic-leading stimulation
    %------------------------------
    fprintf(' | cathodic \n')

    % Determine filenames according to experiment to load and load raw data
    for el = 1:16
        if  strcmp(cfg.Anoexps{exps},'GP09F16CL') || strcmp(cfg.Anoexps{exps},'GP13I16CL')
            file = sprintf('%sctx_pos_%d_layering_100%02d_0001.mat',pth,position,cfg.channelidx(el)-1);
        elseif strcmp(cfg.Anoexps{exps},'GP14F16CL')
            file = sprintf('%sctx_pos%d_layering_100%02d_0001.mat',pth,position,cfg.channelidx(el)-1);
        elseif strcmp(cfg.Anoexps{exps},'GP11K15CL')
            file = sprintf('%sctx_pos%d_layering_mono_100%02d0001.mat',pth,position,tmpchannelidx(el)-1);
        else
            file = sprintf('%sctx_pos%d_layering_100%02d_32_0001.mat',pth,position,cfg.channelidx(el)-1);
        end
        data(el,:,:,:) = Load_AO(file,1,'sweepBefore',50,'sweepAfter',200,'sweeps',32);
    end
    
    % Stimulation artefact blanking
    for el = 1:16
        dataBlanked(el,:,:,:) = ...
            blankstimulus(squeeze(data(el,:,:,:)),22000,50,53);
    end
    clear data
    
    % Filter for MUA
    for el = 1:size(dataBlanked,1)
        for ch = 1:size(dataBlanked,2)
            for sw = 1:size(dataBlanked,3)
                data(el,ch,sw,:) = filtfilt(cfg.filter.bS,cfg.filter.aS,double(squeeze(dataBlanked(el,ch,sw,:))));
            end
        end
    end
    clear el ch sw dataBlanked
    
    % Determine spikes
    for el = 1:16
        [~,MUASpikesC(exps,el,:,:,:)] = spikes(squeeze(data(el,:,:,:)),4);
    end
    MUAC(exps,:,:,:,:) = data;
    clear data

end

% Re-order channels according to Neuronexus electrode mapping
for e = 1:8
    if strcmp(cfg.Anoexps{e},'GP11K15CL')
        MUASpikesA(e,:,:,:,:) = MUASpikesA(e,:,tmpchannelidx,:,:);
        MUASpikesC(e,:,:,:,:) = MUASpikesC(e,:,tmpchannelidx,:,:);
        MUAA(e,:,:,:,:) = MUAA(e,:,tmpchannelidx,:,:);
        MUAC(e,:,:,:,:) = MUAC(e,:,tmpchannelidx,:,:);
    else
        MUASpikesA(e,:,:,:,:) = MUASpikesA(e,:,cfg.channelidx(1:16),:,:);
        MUASpikesC(e,:,:,:,:) = MUASpikesC(e,:,cfg.channelidx(1:16),:,:);
        MUAA(e,:,:,:,:) = MUAA(e,:,cfg.channelidx(1:16),:,:);
        MUAC(e,:,:,:,:) = MUAC(e,:,cfg.channelidx(1:16),:,:);
    end
end
clear tmpchannelidx