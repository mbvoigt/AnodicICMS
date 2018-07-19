%% Load CSD acoustic
%
% Subfunction for Anodic_Cathodic_main
%
% Loads current-source density data for acoustically stimulated trials
%
function csd = Load_CSD_Acoustic(cfg)

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
    
    % Calculate current-source densities
    csdtemp = zeros(size(data,1),size(data,2),size(data,3));
    for el = 1:size(data,1)
        csdtemp(el,:,:) = calculate_csd(squeeze(mean(data(el,:,:,:),3)),0.15,0,[]);
    end
    clear el data
    
    % Baseline correction
    csd = zeros(size(csdtemp,1),size(csdtemp,2),size(csdtemp,3));
    for el = 1:size(csdtemp,1)
        for ch = 1:size(csdtemp,2)
            csd(exps,el,ch,:) = csdtemp(el,ch,:)-mean(csdtemp(el,ch,cfg.BaselineWindow));
        end
    end
    clear el ch csdtemp
    
end