%% Load CSD electric - electrode 9
%
% Subfunction for Anodic_Cathodic_main
%
% Loads current-source density data for electrically stimulated trials
%
% - Stimulation on electrode 9 -
%
function [csdCath, csdAnod, StimCurrent_Cathodic, StimCurrent_Anodic] = Load_CSD_Electric_IO_Electrode_9(cfg)

% preallocation
csdAnod = zeros(length(cfg.SSAnoexps),15,16,5501);
StimCurrent_Anodic = zeros(length(cfg.SSAnoexps),15,30,2201);
csdCath = zeros(length(cfg.SSAnoexps),15,16,5501);
StimCurrent_Cathodic = zeros(length(cfg.SSAnoexps),15,30,2201);

% repeat for every experiment
for exps = 1:length(cfg.SSAnoexps)
    
    % Define folder where data can be found for this specific experiment
    fprintf('%2d | %2d - %s | El. 9 | anodic',exps,length(cfg.SSAnoexps),cfg.SSAnoexps{exps})
    pth = ['K:\ICMS-GP\Data\' cfg.SSAnoexps{exps} '\AlphaOmega\MATFiles\'];
    position = 1; % if more than one position would have been measured this could be implemented here
    
    %----------------------------
    % Anodic-leading stimulation
    %----------------------------
    
    % Determine filenames according to experiment to load and load raw data
    if strcmp(cfg.SSAnoexps{exps},'GP29E15CL') || strcmp(cfg.SSAnoexps{exps},'dummy')
        file = {sprintf('%sctx_pos%d_ss_inf_anodic0001.mat',pth,position), sprintf('%sctx_pos%d_ss_inf_anodic0002.mat',pth,position)};
        data(:,:,:,:) = Load_AO(file,15,'sweepBefore',50,'sweepAfter',200,'sweeps',32);
    else
        file = sprintf('%sctx_pos%d_ss_inf_anodic0001.mat',pth,position);
        data(:,:,:,:) = Load_AO(file,15,'sweepBefore',50,'sweepAfter',200,'sweeps',32);
    end
    
    % Stimulation artefact blanking
    dataBlanked = zeros(size(data,1),size(data,2),size(data,3),size(data,4));
    for str = 1:15
        dataBlanked(str,:,:,:) = ...
            blankstimulus(squeeze(data(str,:,:,:)),22000,50,53);
    end
    clear data
    
    %{
    % filter for MUA
    for el = 1:size(dataBlanked,1)
        for ch = 1:size(dataBlanked,2)
            for sw = 1:size(dataBlanked,3)
                data(el,ch,sw,:) = filtfilt(cfg.filter.bS,cfg.filter.aS,double(squeeze(dataBlanked(el,ch,sw,:))));
            end
        end
    end
    clear el ch sw dataBlanked
    
    % Spikes
    for el = 1:15
        MUA(exps,el,:,:,:) = spikes(squeeze(data(el,:,:,:)));
    end
    clear data
    %}
    
    % Filter for LFP
    for el = 1:size(dataBlanked,1)
        for ch = 1:size(dataBlanked,2)
            for sw = 1:size(dataBlanked,3)
                data(el,ch,sw,:) = filtfilt(cfg.filter.b,cfg.filter.a,double(squeeze(dataBlanked(el,ch,sw,:))));
            end
        end
    end
    clear el ch sw dataBlanked
    
    % Calculate current-source densities
    csd = zeros(size(data,1),size(data,2),size(data,4));
    for el = 1:size(data,1)
        csd(el,:,:) = calculate_csd(squeeze(mean(data(el,:,:,:),3)),0.15,1,9);
    end
    clear el data
    
    % Baseline correction
    for el = 1:size(csd,1)
        for ch = 1:size(csd,2)
            csdAnod(exps,el,ch,:) = csd(el,ch,:)-mean(csd(el,ch,cfg.BaselineWindow));
        end
    end
    clear el ch csd
    
    % Load stimulation current data
    if strcmp(cfg.SSAnoexps{exps},'GP29E15CL') || strcmp(cfg.SSAnoexps{exps},'dummy')
        file = {sprintf('%sctx_pos%d_ss_inf_anodic0001.mat',pth,position), sprintf('%sctx_pos%d_ss_inf_anodic0002.mat',pth,position)};
        dataCurr(:,:,:,:) = Load_AO(file,15,'sweepBefore',50,'sweepAfter',50,'sweeps',32,'probe',NaN);
    else
        file = sprintf('%sctx_pos%d_ss_inf_anodic0001.mat',pth,position);
        dataCurr(:,:,:,:) = Load_AO(file,15,'sweepBefore',50,'sweepAfter',50,'sweeps',32,'probe',NaN);
    end
    StimCurrent_Anodic(exps,:,:,:) = squeeze(double(dataCurr(:,1,:,:))/100);
    clear dataCurr
    
    %------------------------------
    % Cathodic-leading stimulation
    %------------------------------
    fprintf(' | cathodic \n')
    
    % Determine filenames according to experiment to load and load raw data
    if strcmp(cfg.SSAnoexps{exps},'GP29E15CL') || strcmp(cfg.SSAnoexps{exps},'dummy')
        file = {sprintf('%sctx_pos%d_ss_inf0001.mat',pth,position), sprintf('%sctx_pos%d_ss_inf0002.mat',pth,position)};
        data(:,:,:,:) = Load_AO(file,15,'sweepBefore',50,'sweepAfter',200,'sweeps',32);
    else
        file = sprintf('%sctx_pos%d_ss_inf0001.mat',pth,position);
        data(:,:,:,:) = Load_AO(file,15,'sweepBefore',50,'sweepAfter',200,'sweeps',32);
    end
    
    % Stimulation artefact blanking
    for str = 1:15
        dataBlanked(str,:,:,:) = ...
            blankstimulus(squeeze(data(str,:,:,:)),22000,50,53);
    end
    clear data
    
    %{
    % filter for MUA
    for el = 1:size(dataBlanked,1)
        for ch = 1:size(dataBlanked,2)
            for sw = 1:size(dataBlanked,3)
                data(el,ch,sw,:) = filtfilt(cfg.filter.bS,cfg.filter.aS,double(squeeze(dataBlanked(el,ch,sw,:))));
            end
        end
    end
    clear el ch sw dataBlanked
    
    % Spikes
    for el = 1:15
        MUA(exps,el,:,:,:) = spikes(squeeze(data(el,:,:,:)));
    end
    clear data
    %}
    
    % Filter for LFP
    for el = 1:size(dataBlanked,1)
        for ch = 1:size(dataBlanked,2)
            for sw = 1:size(dataBlanked,3)
                data(el,ch,sw,:) = filtfilt(cfg.filter.b,cfg.filter.a,double(squeeze(dataBlanked(el,ch,sw,:))));
            end
        end
    end
    clear el ch sw dataBlanked
    
    % Calculate current-source densities
    for el = 1:size(data,1)
        csd(el,:,:) = p_calculate_csd(squeeze(mean(data(el,:,:,:),3)),0.15,1,9);
    end
    clear el data
    
    % Baseline correction
    for el = 1:size(csd,1)
        for ch = 1:size(csd,2)
            csdCath(exps,el,ch,:) = csd(el,ch,:)-mean(csd(el,ch,cfg.BaselineWindow));
        end
    end
    clear el ch csd
    
    % Load stimulation current data
    if strcmp(cfg.SSAnoexps{exps},'GP29E15CL') || strcmp(cfg.SSAnoexps{exps},'dummy')
        file = {sprintf('%sctx_pos%d_ss_inf0001.mat',pth,position), sprintf('%sctx_pos%d_ss_inf0002.mat',pth,position)};
        dataCurr(:,:,:,:) = Load_AO(file,15,'sweepBefore',50,'sweepAfter',50,'sweeps',32,'probe',NaN);
    else
        file = sprintf('%sctx_pos%d_ss_inf0001.mat',pth,position);
        dataCurr(:,:,:,:) = Load_AO(file,15,'sweepBefore',50,'sweepAfter',50,'sweeps',32,'probe',NaN);
    end
    StimCurrent_Cathodic(exps,:,:,:) = squeeze(double(dataCurr(:,1,:,:))/100);
    
end