%% Load LFP acoustic
%
% Subfunction for Anodic_Cathodic_main
%
% Loads local field potential data for acoustically stimulated trials
%
function lfp = Load_LFP_Acoustic(cfg)

% Variable preallocation
lfp = zeros(length(cfg.SSAnoexps),17,16,30,5501);

% repeat for every experiment
for exps = 1:length(cfg.SSAnoexps)
    
    % Define folder where data can be found for this specific experiment
    fprintf('%2d | %d - %s | acoustic\n',exps,length(cfg.SSAnoexps),cfg.SSAnoexps{exps})
    pth = ['K:\ICMS-GP\Data\' cfg.SSAnoexps{exps} '\AlphaOmega\MATFiles\'];
    position = 1; % if more than one position would have been measured this could be implemented here
    
    % Determine filenames according to experiment to load and load raw data
    if strcmp(cfg.SSAnoexps{exps},'GP01L15CL')
        file = sprintf('%sctx_pos%d_preexp_click0001.mat',pth,position);
        data(:,:,:,:) = Load_AO(file,17,'sweepBefore',50,'sweepAfter',200,'sweeps',30);
    elseif strcmp(cfg.SSAnoexps{exps},'GP29E15CL')
        file = {sprintf('%sctx_pos%d_click0002.mat',pth,position),sprintf('%sctx_pos%d_click0003.mat',pth,position)};
        data(:,:,:,:) = Load_AO(file,17,'sweepBefore',50,'sweepAfter',200,'sweeps',30);
    else
        file = sprintf('%sctx_pos%d_click0001.mat',pth,position);
        data(:,:,:,:) = Load_AO(file,17,'sweepBefore',50,'sweepAfter',200,'sweeps',30);
    end
    
    % Convert raw data to double precision (AlphaOmega default: 16-bit int)
    data = double(data);
    
    % Filter for LFP
    for el = 1:size(data,1)
        for ch = 1:size(data,2)
            for sw = 1:size(data,3)
                data(el,ch,sw,:) = filtfilt(cfg.filter.b,cfg.filter.a,squeeze(data(el,ch,sw,:)));
            end
        end
    end
    clear el ch sw
    
    lfp(exps,:,:,:,:) = data;
    clear data
    
end