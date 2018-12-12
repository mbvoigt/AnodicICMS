%% Differential cortical activation according to the leading-phase polarity
%  of intracortical microstimulation
%
% The present script, inlcuding all accompanying sub-scripts, was used to
% analyze the data and prepare the figures for the manuscript entitled
% "Cathodic-leading pulses are more effective than anodic-leading pulses in
% intracortical microstimulation of the auditory cortex" 
% by Mathias B. Voigt and Andrej Kral.
%
%
% Abbreviations:
%   CSD - Current-source density
%   IO  - Input-Output function
%   LFP - Local field potential
%   MUA - Multi unit activity
%
% External dependencies: 
%   - Matlab                                         (Version used: 9.3)
%   - Signal Processing Toolbox                      (Verison used: 7.5)
%   - Symbolic Math Toolbox                          (Version used: 8.0)
%   - Statistics and Machine Learning Toolboxes      (Version used: 11.2)
%
%
% This work is licensed under the MIT license.
% See the accompanying LICENSE.txt for detail.
%
% Copyright (c) 2018 Mathias Voigt <voigt.mathias@mh-hannover.de>
% Institute of AudioNeuroTechnology, Hannover Medical School
%
%
%% Anodic_Cathodic_main

% Clean up
clear variables
close all
clc

% Logging of script progress
cfg.logfile = ['D:\mbv\log\Manuscript_AnoBip_' datestr(now,'yy-mm-dd_HH_MM_SS') '.txt'];
mbv_log('Preparing Analyses for Manuscript - Anodic vs. Cathodic leading & Monopolar vs. Bipolar','file',cfg.logfile);

% Find experiment IDs which were stimulated both anodic- and
% cathodic-leading
cfg.Anoexps = findexperiments('expression','*layer*anod*');
mbv_log(sprintf('Number of anodic experiments to analyze: %2d',length(cfg.Anoexps)),'file',cfg.logfile);
mbv_log(['Anodic experiments to analyze: ' strjoin(cfg.Anoexps,', ')],'file',cfg.logfile);

% Find experiment IDs which were stimulated anodic-leading with varying current
cfg.SSAnoexps = findexperiments('expression','*SS*anod*');
mbv_log(sprintf('Number of input/output experiments to analyze: %2d',length(cfg.SSAnoexps)),'file',cfg.logfile);
mbv_log(['Input/output experiments to analyze: ' strjoin(cfg.SSAnoexps,', ')],'file',cfg.logfile);

% Number of channels to be analyzed (1st 16 channels of the double shank
% electrodes)
cfg.noChannels = 16;
mbv_log(['Number of channels: ' num2str(cfg.noChannels)],'file',cfg.logfile);

%
% Channel mapping from the Neuronexus electrode array contacts to the 
% Alpha Omega recording channels
%
cfg.channelidx = [9,8,10,7,11,6,12,5,13,4,14,3,15,2,16];
%1,25,24,26,23,27,22,28,21,29,20,30,19,31,18,32,17,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48];
mbv_log(['Channel order used: ' sprintf('%02d, ',cfg.channelidx)],'file',cfg.logfile);

% Time windows for linear interpolation to blank the electrical stimulation
% artefact, 'blf' blanking start time [ms], 'blt' blanking stop time [ms]
cfg.blf = 50;
cfg.blt = 53;
mbv_log(['Blanking interval starts: ' sprintf('%02d ms',cfg.blf)],'file',cfg.logfile);
mbv_log(['Blanking interval stops: ' sprintf('%02d ms',cfg.blt)],'file',cfg.logfile);

% Time window for epoch cutting of continuous recordings
% Here: 50 ms before and 200 ms after stimulation trigger timestamp
cfg.sweepBefore = 50;
cfg.sweepAfter = 200;
mbv_log(sprintf('Analyzing %3d ms before and %3d ms after trigger',cfg.sweepBefore,cfg.sweepAfter),'file',cfg.logfile);

% Baseline time window
cfg.BaselineWindow = 1:cfg.sweepBefore*22;
mbv_log(sprintf('Baseline window: 0 ms : %3d ms',cfg.BaselineWindow(end)/22),'file',cfg.logfile);

%
% Filter settings:
%
% LFP
[z,p,k] = butter(2,150*2/22000,'low');
[cfg.filter.b,cfg.filter.a] = zp2tf(z,p,k);
clear z p k
% MUA
[z,p,k] = butter(2,300*2/22000,'high');
[cfg.filter.bS,cfg.filter.aS] = zp2tf(z,p,k);
clear z p k

%% Acoustic stimulation - Calibration of acoustic click intensity

% Load calibration data - 50 µs click, 120:-5:40 dB_Attenuation
data = Load_AudiologyLab('K:\ICMS-GP\Data\Calibration_2016_01_29\click_120-40dBAtt.001',1);
data_calib = squeeze(data.data(:,1,:,:));
data_calib = permute(data_calib,[3 2 1]);

% peak-2-peak values for different attenuations of 50 µs clicks
calib_click = range(mean(data_calib,2),3);

% Load calibration data - 15 ms tone, 0:5:80 dB_Attenuation
data = Load_AudiologyLab('K:\ICMS-GP\Data\Calibration_2016_01_29\tone15ms_0-80dBSPL.001',1);
data_calib = squeeze(data.data(:,1,:,:));
data_calib = permute(data_calib,[3 2 1]);

% peak-2-peak values for different attenuations of 15 ms tone bursts
calib_tone = range(mean(data_calib,2),3);

click_att=120:-5:40;
fprintf('%2d dB attenuation equals 70 dB SPL in peak-to-peak amplitude\n',click_att((find(calib_tone(15) >= calib_click,1,'last')+1)));

dBAtteqto70dBSPL = click_att((find(calib_tone(15) >= calib_click,1,'last')+1));
lowestAttenuation = 70-(120-dBAtteqto70dBSPL);
highestAttenuation = 70+(dBAtteqto70dBSPL-40);

% Acoustic click intensity in dB SPL_peakequivalent
acoustic_strength = lowestAttenuation:5:highestAttenuation; % 15:5:95 dB SPL_pe
fprintf('Acoustic click intensity: %s dB SPL_pe\n',sprintf('%d ',acoustic_strength))
clear data data_calib calib_click calib_tone click_att dBAtteqto70dBSPL lowestAttenuation highestAttenuation

%% Acoustic stimulation - Load acoustic click data

% Load acoustically stimulated LFP data
lfp = Load_LFP_Acoustic(cfg);

% Figure 2A - LFP profile
fH = figure();
hold on
for ch = 1:16
    plot(squeeze(mean(lfp(:,ch,:)))-(200*ch),'k')
end
plot([1100 1100],[-17*4000 0],'--k')
set(gca,'YLim',[-17*200 0],'XLim',[880 2201])
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\Acoustic_LFP_example_profile')
close(fH)
clear lfp

% Load acoustically stimulated CSD data
csd = Load_CSD_Acoustic(cfg);

% Figure 2A - CSD profile
fH = figure();
hold on
for ch = 1:16
    plot(squeeze(mean(csd(:,ch,:)))-(4000*ch),'k')
    fillhelp=squeeze(mean(csd(:,ch,:)));
    fillhelp(fillhelp<0) = 0;
    patch([1,1:5501,5501],[-4000*ch;fillhelp-(4000*ch);-4000*ch],[0 0 0],'FaceColor','k','EdgeColor','none')
    plot([1 5501],[-4000*ch -4000*ch],'--k')
end
plot([1100 1100],[-17*4000 0],'--k')
plot([1000 1000],[-9000 -5000],'k','LineWidth',2)
set(gca,'YLim',[-17*4000 0],'XLim',[880 2201])
set(gcf,'Renderer','painters')
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\Acoustic_CSD_example_profile')
close(fH)

% Determine (sink) peak amplitudes in the time window from 0 to 50 ms post-stimulation
amp_click = squeeze(max(csd(:,:,:,1100:2200),[],4));
clear csd

% Save result of calculations for convenience
save('D:\mbv\temp\Manuscript_anodic\Acoustic_CSDs','amp_click','cfg')
%load('D:\mbv\temp\Manuscript_anodic\Acoustic_CSDs','amp_click')

%% Acoustic stimulation - Input/Output function
% Determine Input-Output functions and non-linear fits of acoustic data
[beta_click, beta_clickM, Rsquared_click] = Analyze_IO_curves_acoustic(cfg,amp_click,0);
% Figures 4C, 4D

%% Electric stimulation - Input/Output function - Load data & start processing

%----------------------------
% Stimulation on electrode 1
%----------------------------

% Load CSD data
[csdCath, csdAnod, StimCurrent_Cathodic, StimCurrent_Anodic] = Load_CSD_Electric_IO_Electrode_1(cfg);

% Save CSDs to file for convenience
% save('D:\mbv\temp\Manuscript_anodic\IO_anocath_El1','csdCath','csdAnod','StimCurrent_Cathodic','StimCurrent_Anodic')
% load('D:\mbv\temp\Manuscript_anodic\IO_anocath_El1','csdCath','csdAnod','StimCurrent_Cathodic','StimCurrent_Anodic')

cfg.saveFileString = 'El1';

% Calibrate the stimulation currents
strengths = StimulationCurrentCalibration(StimCurrent_Cathodic,StimCurrent_Anodic);
% Figures S1A to S1E

% Stimulation strengths saved for convenience:
% strengths = [0.144000000000000 0.286000000000000 0.640888888888889 1.10972222222222 1.70572222222222 2.04505555555556 2.49694444444445 2.88944444444444 3.40527777777778 3.88100000000000 4.43355555555556 5.02438888888889 5.59955555555556 15.0966666666667 36.7231666666667];

% Determine stimulation thresholds
[thrC_E1,thrA_E1] = Analyze_IO_thresholds_single(csdAnod,csdCath,strengths);

% Peak (sink) amplitudes
ampA = squeeze(max(csdAnod(:,:,:,1100:2200),[],4));
ampC = squeeze(max(csdCath(:,:,:,1100:2200),[],4));

% Input/Output function and non-linear fits of electric stimulation
[~,~,RsquaredC_E1,RsquaredA_E1,betaMC,betaMA,~,~] = Analyze_IO_curves(cfg,ampA,ampC);
clear RsquaredC RsquaredA RsquaredMC RsquaredMA StimCurrent_Cathodic StimCurrent_Anodic

% Plot IO figures
Plot_IO_figures(cfg,ampC,ampA,strengths,betaMC,betaMA)
% Figures 4,S3A, S3B

% Determine dynamic range
[DRC_E1,DRA_E1,DRClick] = Analyze_IO_curves_DynamicRange(cfg,strengths,ampA,ampC,amp_click);
% Figure 4E

% Re-save amplitude data to different variables for later comparisons
ampA_E1 = ampA;
ampC_E1 = ampC;
clear ampA ampC

%----------------------------
% Stimulation on electrode 9
%----------------------------

% Load CSD data
[csdCath, csdAnod, ~, ~] = Load_CSD_Electric_IO_Electrode_9(cfg);

% Save CSDs to file for convenience
% save('D:\mbv\temp\Manuscript_anodic\IO_anocath_El9','csdCath','csdAnod','StimCurrent_Cathodic','StimCurrent_Anodic')
% load('D:\mbv\temp\Manuscript_anodic\IO_anocath_El9','csdCath','csdAnod','StimCurrent_Cathodic','StimCurrent_Anodic')

cfg.saveFileString = 'El9';

% Determine stimulation thresholds
[thrC_E9,thrA_E9] = Analyze_IO_thresholds_single(csdAnod,csdCath,strengths);

% Peak (sink) amplitudes
ampA = squeeze(max(csdAnod(:,:,:,1100:2200),[],4));
ampC = squeeze(max(csdCath(:,:,:,1100:2200),[],4));

% Input/Output function and non-linear fits of electric stimulation
[~,~,RsquaredC_E9,RsquaredA_E9,betaMC,betaMA,~,~] = Analyze_IO_curves(cfg,ampA,ampC);
clear RsquaredC RsquaredA RsquaredMC RsquaredMA StimCurrent_Cathodic StimCurrent_Anodic

% Plot IO figures
% Attention: Change YLimit in script for stimulation on electrode 9!
Plot_IO_figures(cfg,ampC,ampA,strengths,betaMC,betaMA)
% Figures S3A, S3B

% Determine dynamic range
[DRC_E9,DRA_E9,~] = Analyze_IO_curves_DynamicRange(cfg,strengths,ampA,ampC,amp_click);
% Figure S4A

% Re-save amplitude data to different variables for later comparisons
ampA_E9 = ampA;
ampC_E9 = ampC;
clear ampA ampC

%-----------------------------
% Stimulation on electrode 16
%-----------------------------

% Load CSD data
[csdCath, csdAnod, ~, ~] = Load_CSD_Electric_IO_Electrode_16(cfg);

% Save CSDs to file for convenience
% save('D:\mbv\temp\Manuscript_anodic\IO_anocath_El16','csdCath','csdAnod','StimCurrent_Cathodic','StimCurrent_Anodic')
% load('D:\mbv\temp\Manuscript_anodic\IO_anocath_El16','csdCath','csdAnod','StimCurrent_Cathodic','StimCurrent_Anodic')

cfg.saveFileString = 'El16';

% Determine stimulation thresholds
[thrC_E16,thrA_E16] = Analyze_IO_thresholds_single(csdAnod,csdCath,strengths);

% Peak (sink) amplitudes
ampA = squeeze(max(csdAnod(:,:,:,1100:2200),[],4));
ampC = squeeze(max(csdCath(:,:,:,1100:2200),[],4));

% Input/Output function and non-linear fits of electric stimulation
[betaC,betaA,RsquaredC_E16,RsquaredA_E16,betaMC,betaMA,~,~] = Analyze_IO_curves(cfg,ampA,ampC);
clear RsquaredC RsquaredA RsquaredMC RsquaredMA StimCurrent_Cathodic StimCurrent_Anodic

% Plot IO figures
Plot_IO_figures(cfg,ampC,ampA,strengths,betaMC,betaMA)
% Figures S3A, S3B

% Determine dynamic range
[DRC_E16,DRA_E16,~] = Analyze_IO_curves_DynamicRange(cfg,strengths,ampA,ampC,amp_click);
% Figure S4B

% Re-save amplitude data to different variables for later comparisons
ampA_E16 = ampA;
ampC_E16 = ampC;
clear ampA ampC

%% Electric stimulation - Input/Output function - Stimulation threshold comparison

% Concatenate threshold values for statistical testing
y = [strengths(min(thrC_E1,[],2)),strengths(min(thrA_E1,[],2)),...
    strengths(min(thrC_E9,[],2)),strengths(min(thrA_E9,[],2)),...
    strengths(min(thrC_E16,[],2)),strengths(min(thrA_E16,[],2))];
% Grouping variables
f1 = [ones(1,16),ones(1,16)*2,ones(1,16)*3]; % Stimulation electrode (1, 9, 16)
f2 = repmat([ones(1,8),ones(1,8)*2],1,3); % Stimulus polarity (cathodic-leading, anodic-leading)
group = {f1,f2};

% 2-way ANOVA
% [p,tbl,stats] = anovan(y,group,'model','interaction','varnames',{'Electrode','Polarity'})

% post-hoc Tukey's HSD
% c = multcompare(stats)

% Figure S2
fH = figure();
boxplot(y,group,'Colors','k')
set(gca,'YLim',[0 6])
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\IO_anocath_threshold_boxplot')
close(fH)

%% Electric stimulation - Input/Output function - Electrode comparisons

% Standard-deviation of Input/Output functions
% Figure S3C
fH = figure();
hold on
plot(strengths,(std(mean(ampC_E1,3),1)),'Color',[0 0 120/255])
plot(strengths,(std(mean(ampA_E1,3),1)),'Color',[120/255 0 0])
plot(strengths,(std(mean(ampC_E9,3),1)),'--','Color',[0 0 120/255])
plot(strengths,(std(mean(ampA_E9,3),1)),'--','Color',[120/255 0 0])
plot(strengths,(std(mean(ampC_E16,3),1)),':','Color',[0 0 120/255])
plot(strengths,(std(mean(ampA_E16,3),1)),':','Color',[120/255 0 0])
set(gca,'Xscale','log')
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\IO_CSD_Anodic_singleExperiments_standarddeviation')
close(fH)

% Goodness of fit - R^2
% Figure S3D
fH = figure();
hold on
plot(ones(1,8)*1,RsquaredC_E1,'o','MarkerEdgeColor','none','MarkerFaceColor',[0 0 120/255])
plot(ones(1,8)*2,RsquaredA_E1,'o','MarkerEdgeColor','none','MarkerFaceColor',[120/255 0 0])
plot(ones(1,8)*3,RsquaredC_E9,'o','MarkerEdgeColor','none','MarkerFaceColor',[0 0 120/255])
plot(ones(1,8)*4,RsquaredA_E9,'o','MarkerEdgeColor','none','MarkerFaceColor',[120/255 0 0])
plot(ones(1,8)*5,RsquaredC_E16,'o','MarkerEdgeColor','none','MarkerFaceColor',[0 0 120/255])
plot(ones(1,8)*6,RsquaredA_E16,'o','MarkerEdgeColor','none','MarkerFaceColor',[120/255 0 0])
plot(ones(1,8)*7,Rsquared_click,'o','MarkerEdgeColor','none','MarkerFaceColor',[0 0 0])
boxplot([RsquaredC_E1' RsquaredA_E1' RsquaredC_E9' RsquaredA_E9' RsquaredC_E16' RsquaredA_E16' Rsquared_click'],'Color','k')
set(gca,'XLim',[0 8],'YLim',[0.7 1])
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\IO_CSD_Rsquared_comparison')
close(fH)

% Statistics
% [p,tbl,stats]=kruskalwallis([RsquaredC_E1' RsquaredA_E1' RsquaredC_E9' RsquaredA_E9' RsquaredC_E16' RsquaredA_E16' Rsquared_click'])
% c = multcompare(stats) % post-hoc Tukey HSD

% Difference Cathodic-Anodic
% Figure 4B
fH = figure();
hold on
scatter(ones(8,1),mean(mean(ampC_E1,3),2)-mean(mean(ampA_E1,3),2),'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.2)
scatter(ones(8,1)*2,mean(mean(ampC_E9,3),2)-mean(mean(ampA_E9,3),2),'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.2)
scatter(ones(8,1)*3,mean(mean(ampC_E16,3),2)-mean(mean(ampA_E16,3),2),'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.2)
boxplot([mean(mean(ampC_E1,3),2)-mean(mean(ampA_E1,3),2),mean(mean(ampC_E9,3),2)-mean(mean(ampA_E9,3),2),mean(mean(ampC_E16,3),2)-mean(mean(ampA_E16,3),2)])
plot([0 4],[0 0],':k')
set(gca,'XLim',[0 4],'YLim',[-1000 2500])
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\IO_CSD_Diff')
close(fH)

% Dynamic range
% Figure 4F
fH = figure();
hold on
plot(ones(1,8),DRC_E1,'ok')
plot(ones(1,8)*2,DRA_E1,'ok')
plot(ones(1,8)*3,DRC_E9,'ok')
plot(ones(1,8)*4,DRA_E9,'ok')
plot(ones(1,8)*5,DRC_E16,'ok')
plot(ones(1,8)*6,DRA_E16,'ok')
plot(ones(1,8)*7,DRClick,'ok')
boxplot([DRC_E1' DRA_E1' DRC_E9' DRA_E9' DRC_E16' DRA_E16' DRClick'],'Colors','k');
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\IO_CSD_DynamicRange')
close(fH)
% Statistics
% [p,tbl,stats] = kruskalwallis([DRC_E1' DRA_E1' DRC_E9' DRA_E9' DRC_E16' DRA_E16' DRClick']);
% c = multcompare(stats);

fprintf('DR - cathodic: %5.3f +- %5.3f\n',mean(nanmean([DRC_E1' DRC_E9' DRC_E16'])),std(nanmean([DRC_E1' DRC_E9' DRC_E16'])))
fprintf('DR - cathodic: %5.3f +- %5.3f\n',mean(nanmean([DRA_E1' DRA_E9' DRA_E16'])),std(nanmean([DRA_E1' DRA_E9' DRA_E16'])))
fprintf('DR - click:    %5.3f +- %5.3f\n',nanmean(DRClick),nanstd(DRClick))

%% Electric stimulation - Input/Output function - Statistics - 3-way ANOVA

% Take mean over all recording electrodes
ampA_E1 = squeeze(mean(ampA_E1,3));
ampC_E1 = squeeze(mean(ampC_E1,3));
ampA_E9 = squeeze(mean(ampA_E9,3));
ampC_E9 = squeeze(mean(ampC_E9,3));
ampA_E16 = squeeze(mean(ampA_E16,3));
ampC_E16 = squeeze(mean(ampC_E16,3));

% Concatenate for testing
values = [ampC_E1(:);ampC_E9(:);ampC_E16(:);ampA_E1(:);ampA_E9(:);ampA_E16(:)];

% Define factors
gStr = repmat((kron((1:15)',ones(8,1))),6,1); % Stimulation strength
gEle = repmat((kron((1:3)',ones(8*15,1))),2,1); % Stimulated electrode
gPol = [ones(8*15*3,1);ones(8*15*3,1)*2];  % Leading-phase polarity

% ANOVA
[p,tbl,~]=anovan(values,{gStr,gEle,gPol},'varnames',{'strength','electrode','polarity'},'model','full');

% Calculate R^2 as ratio of variance explained by the different factors
SSstr = tbl{2,2}/tbl{10,2};
SSEle = tbl{3,2}/tbl{10,2};
SSPol = tbl{4,2}/tbl{10,2};

fprintf('3-way ANOVA:\n p = %6.4f\n R^2 strength: %d\n R^2 electrode: %d\n R^2 polarity: %d\n',p,SSstr,SSEle,SSPol);

%% Electric stimulation - Varying depth - LFP - Load data

% Load current-source density data for stimulation with varying depth
[lfpAnod,lfpCath] = Load_LFP_Electric_VarDepth(cfg);

% Save LFP data to file for convenience
% save('D:\mbv\temp\Manuscript_anodic\anocath_layering_LFPs','lfpAnod','lfpCath','-v7.3')
% load('D:\mbv\temp\Manuscript_anodic\anocath_layering_LFPs','lfpAnod','lfpCath')

%% Electric stimulation - Varying depth - LFP - Figures

% Average over stimulus repetitions for a stimulation current of ~ 6 µA
lfpC = squeeze(mean(lfpCath(:,9,:,:,:),4));
lfpA = squeeze(mean(lfpAnod(:,9,:,:,:),4));
clear lfpCath lfpAnod

% Figure 2B - Cathodic
fH = figure();
hold on
for ch = 1:16
    plot(squeeze(mean(lfpC(:,ch,:),1))-(200*ch),'k')
end
plot([1100 1100],[-17*4000 0],'--k')
set(gca,'YLim',[-17*200 0],'XLim',[880 2201])
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\Electric_LFP_example_profile_cath')
close(fH)

% Figure 2B - Anodic
fH = figure();
hold on
for ch = 1:16
    plot(squeeze(mean(lfpA(:,ch,:),1))-(200*ch),'k')
end
plot([1100 1100],[-17*4000 0],'--k')
set(gca,'YLim',[-17*200 0],'XLim',[880 2201])
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\Electric_LFP_example_profile_anod')
close(fH)

%
% Figure 3A - "Activity matrices"
% 
ampAnod = squeeze(range(lfpA(:,:,:,1100:2200),4));
ampCath = squeeze(range(lfpC(:,:,:,1100:2200),4));
clear lfpA lfpC

% Colorbar
fH = figure();
caxis([-600 600])
colormap(bluewhitered())
colorbar()
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\AM_LFP_colorbar')
close(fH)

% Anodic
fH = figure();
imagesc(squeeze(mean(ampAnod))')
caxis([-600 600])
colormap(bluewhitered())
axis square
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\AM_LFP_Anodic')
close(fH)

% Cathodic
fH = figure();
imagesc(squeeze(mean(ampCath))')
caxis([-600 600])
colormap(flip(bluewhitered()))
axis square
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\AM_LFP_Cathodic')
close(fH)

% Difference
fH = figure();
imagesc((squeeze(mean(ampCath))-squeeze(mean(ampAnod)))')
caxis([-600 600])
colormap(flip(bluewhitered()))
axis square
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\AM_LFP_Diff')
close(fH)

clear ampAnod ampCath

%% Electric stimulation - Varying depth - CSD - Load data

% Load current-source density data for stimulation with varying depth
[csdAnod,csdCath] = Load_CSD_Electric_VarDepth(cfg);

% Save CSD data to file for convenience
% save('D:\mbv\temp\Manuscript_anodic\anocath_layering_CSDs','csdAnod','csdCath')
% load('D:\mbv\temp\Manuscript_anodic\anocath_layering_CSDs','csdAnod','csdCath')

%% Electric stimulation - Varying depth - CSD - Example pictures

% CSD examples
% Figure 2C, 2D, 2E
% Examples used: Electrode 1, 9, and 16
for el = 1:16
    fH = figure();
    hold on
    for ch = 1:16
        plot(squeeze(mean(csdCath(:,el,ch,880:2200)))-(4000*ch),'k')
        plot([1 1321],[-4000*ch -4000*ch],'k:')
    end
    plot([220 220],[-17*4000 0],':k')
    plot([10 10],[-5000 0],'k','LineWidth',2)
    set(gca,'XLim',[0 1321],'XTick',0:220:1500,'XTickLabel',-10:10:200)
    print(fH,'-dpdf','-r1200',['D:\mbv\temp\Manuscript_anodic\plots\CSD\CSD_Ano_CSD_GrandMean_Cath_stimEl' num2str(el)])
    close(fH)
    
    fH = figure();
    hold on
    for ch = 1:16
        plot(squeeze(mean(csdAnod(:,el,ch,880:2200)))-(4000*ch),'k')
        plot([1 1321],[-4000*ch -4000*ch],'k:')
    end
    plot([220 220],[-17*4000 0],':k')
    plot([10 10],[-5000 0],'k','LineWidth',2)
    set(gca,'XLim',[0 1321],'XTick',0:220:1500,'XTickLabel',-10:10:200)
    print(fH,'-dpdf','-r1200',['D:\mbv\temp\Manuscript_anodic\plots\CSD\CSD_Ano_CSD_GrandMean_Anod_stimEl' num2str(el)])
    close(fH)
end

%% Electric stimulation - Varying depth - CSD - Sink analysis

% We analyzed only sinks, therefore all negative values were set to 0:
sinkCath = csdCath;
sinkCath(sinkCath<0)=0;
sinkAnod = csdAnod;
sinkAnod(sinkAnod<0)=0;

% Determination of peak amplitudes
ampAnod = max(sinkAnod(:,:,:,1100:2200),[],4);
ampCath = max(sinkCath(:,:,:,1100:2200),[],4);

% Figure 3B - "Activity matrices"

% Colorbar
fH = figure();
caxis([-35000 35000])
colormap(bluewhitered())
colorbar()
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\AM_colorbar')
close(fH)

% Anodic
fH = figure();
imagesc(squeeze(mean(ampAnod))')
caxis([-30000 30000])
colormap(bluewhitered())
axis square
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\AM_Anod')
close(fH)

% Cathodic
fH = figure();
imagesc(squeeze(mean(ampCath))')
caxis([-30000 30000])
colormap(flip(bluewhitered()))
axis square
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\AM_Cath')
close(fH)

% Difference
fH = figure();
imagesc((squeeze(mean(ampCath))-squeeze(mean(ampAnod)))')
caxis([-35000 35000])
colormap(flip(bluewhitered()))
axis square
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\AM_Diff')
close(fH)

% Figure 3C
fH = figure();
hold on
errorbar(mean(mean(ampAnod,3)),16:-1:1,std(mean(ampAnod,3))./sqrt(8),'o-r','horizontal')
errorbar(mean(mean(ampCath,3)),16:-1:1,std(mean(ampCath,3))./sqrt(8),'o-k','horizontal')
errorbar(mean(mean(mean(ampAnod,3))),0,std(mean(mean(ampAnod,3)))./sqrt(8),'o-r','horizontal')
errorbar(mean(mean(mean(ampCath,3))),-1,std(mean(mean(ampCath,3)))./sqrt(8),'o-k','horizontal')
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\CSD_Ano_CSDamplitudes')
close(fH)

% Measures:
fprintf('Cathodic, mean: %3.2f +- %3.2f\n',mean(mean(mean(ampCath,3)))/1000,std(mean(mean(ampCath,3)))/1000)
fprintf('Anodic, mean: %3.2f +- %3.2f\n',mean(mean(mean(ampAnod,3)))/1000,std(mean(mean(ampAnod,3)))/1000)

% Statistic
% Difference in grand mean
% [p,h,stats]=signrank(mean(mean(ampAnod,3)),mean(mean(ampCath,3)))

% repeated measures ANOVA
anim = categorical(repmat((1:8)',16,1));
stimel = categorical(reshape(repmat(1:16,8,1),[],1));
t = table(anim,stimel,reshape(mean(ampCath,3),[],1),reshape(mean(ampAnod,3),[],1),'VariableNames',{'animal','StimEl','dC','dA'});
within = table(categorical([1 2]',[1 2],{'cathodic','anodic'}),'VariableNames',{'polarity'});
rm = fitrm(t,'dC,dA ~ StimEl','WithinDesign',within);
ranovatbl = ranova(rm,'WithinModel','polarity'); %#ok! values are manually examined and reported in the manuscript
c = multcompare(rm,'polarity','By','StimEl'); %#ok! post-hoc Tukey's HSD

%
% PEAK LATENCY
%
% Determination of peak latencies
latAnod = zeros(8,16,16);
latCath = zeros(8,16,16);
for exps = 1:8
    for el = 1:16
        for ch = 1:16
            latAnod(exps,el,ch) = find(sinkAnod(exps,el,ch,1100:2200)==ampAnod(exps,el,ch),1,'first')/22;
            latCath(exps,el,ch) = find(sinkCath(exps,el,ch,1100:2200)==ampCath(exps,el,ch),1,'first')/22;
        end
    end
end

% Figure 3D
fH = figure();
hold on
errorbar(mean(mean(latAnod,3)),16:-1:1,std(mean(latAnod,3))./sqrt(8),'o-r','horizontal')
errorbar(mean(mean(latCath,3)),16:-1:1,std(mean(latCath,3))./sqrt(8),'o-k','horizontal')
errorbar(mean(mean(mean(latAnod,3))),0,std(mean(mean(latAnod,3)))./sqrt(8),'o-r','horizontal')
errorbar(mean(mean(mean(latCath,3))),-1,std(mean(mean(latCath,3)))./sqrt(8),'o-k','horizontal')
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\CSD_Ano_CSDlatencies')
close(fH)

% Measures:
fprintf('Cathodic, mean: %5.2f +- %6.2f\n',mean(mean(mean(latCath,3))),std(mean(mean(latCath,3))))
fprintf('Anodic, mean: %5.2f +- %6.2f\n',mean(mean(mean(latAnod,3))),std(mean(mean(latAnod,3))))

% Statistic
% Difference in grand mean
[p,~,stats]=signrank(mean(mean(latAnod,3)),mean(mean(latCath,3))); %#ok variables need to be inspected manually

% repeated measures ANOVA
anim = categorical(repmat((1:8)',16,1));
stimel = categorical(reshape(repmat(1:16,8,1),[],1));
t = table(anim,stimel,reshape(mean(latCath,3),[],1),reshape(mean(latAnod,3),[],1),'VariableNames',{'animal','StimEl','lC','lA'});
within = table(categorical([1 2]',[1 2],{'cathodic','anodic'}),'VariableNames',{'polarity'});
rm = fitrm(t,'lC,lA ~ StimEl','WithinDesign',within);
ranovatbl = ranova(rm,'WithinModel','polarity'); %#ok! Values are manually examined and reported in the manuscript
c = multcompare(rm,'polarity','By','StimEl'); %#ok! Values are manually examined and reported in the manuscript

% Correlation: Amplitude - Latency
mampAnod = mean(mean(ampAnod,3));
mampCath = mean(mean(ampCath,3));
mlatAnod = mean(mean(latAnod,3));
mlatCath = mean(mean(latCath,3));

[rho,pval]=corr(mampCath(:),mlatCath(:),'type','spearman'); %#ok! Values are manually examined and reported in the manuscript
[rho,pval]=corr(mampAnod(:),mlatAnod(:),'type','spearman');

%% Electric stimulation - Varying depth - CSD - Normalization

% Normalization
sinkCathN = zeros(8,16,16,5501);
sinkAnodN = zeros(8,16,16,5501);
for exps = 1:8
    for el = 1:16
        for ch = 1:16
            sinkCathN(exps,el,ch,:) = sinkCath(exps,el,ch,:)./max(max(sinkCath(exps,el,:,:)));
            sinkAnodN(exps,el,ch,:) = sinkAnod(exps,el,ch,:)./max(max(sinkAnod(exps,el,:,:)));
        end
    end
end

% Normalized CSD sinks example
% Figure 5A
fH = figure();
hold on
for ch = 1:16
    plot(squeeze(sinkCathN(6,9,ch,880:2200))-(1*ch),'k')
%     plot(squeeze(sinkAnodN(6,9,ch,880:2200))-(1*ch),'r')
end
plot([0 1100],[-7.2 -7.2],':k')
plot([220 220],[-17 0],':k')
plot([220 440],[-1.5 -1.5],'k','LineWidth',2) % time axis scale bar
plot([1100 1100],[-2 -1],'k','LineWidth',2) % amplitude axis scale bar
set(gca,'XLim',[1 1320],'YLim',[-17*1 0],'XTick',0:220:1100,'XTickLabel',-10:10:100)
print(fH,'-dsvg','-r1200',['D:\mbv\temp\Manuscript_anodic\plots\CSD_Ano_CSDSinksnorm_' cfg.Anoexps{6} '_example'])
close(fH)

%% Electric stimulation - Varying depth - CSD - Activity matrices

% Threshold value
thr = 0.8;

% Thresholding
sinkCathNthr = any(sinkCathN(:,:,:,1100:2200)>thr,4);
sinkAnodNthr = any(sinkAnodN(:,:,:,1100:2200)>thr,4);

% Difference matrix
dif = double(sinkAnodNthr) - double(sinkCathNthr);

% Plot single experiment data
% Figure 5B, S5
for exps = 1:8
    % Cathodic-leading
    fH = figure();
    imagesc(squeeze(sinkCathNthr(exps,:,:))')
    axis square
    colormap(flipud(gray()))
    print(fH,'-dsvg','-r1200',['D:\mbv\temp\Manuscript_anodic\plots\CSD_Ano_SinkMatrix_' cfg.Anoexps{exps} '_cathodic'])
    close(fH)
    % Anodic-leading
    fH = figure();
    imagesc(squeeze(sinkAnodNthr(exps,:,:))')
    axis square
    colormap(flipud(gray()))
    print(fH,'-dsvg','-r1200',['D:\mbv\temp\Manuscript_anodic\plots\CSD_Ano_SinkMatrix_' cfg.Anoexps{exps} '_anodic'])
    close(fH)
    % Difference
    fH = figure();
    imagesc(squeeze(dif(exps,:,:))')
    axis square
    colormap(bluewhitered()) % The function bluewhitered
    print(fH,'-dsvg','-r1200',['D:\mbv\temp\Manuscript_anodic\plots\CSD_Ano_SinkMatrix_' cfg.Anoexps{exps} '_dif'])
    close(fH)
end

% Plot grand mean activity matrices
% Figure 5C
% Cathodic-leading
fH = figure();
imagesc(squeeze(mean(sinkCathNthr(:,:,:)))')
axis square
colormap(flipud(gray()))
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\CSD_Ano_SinkMatrix_GrandMean_cathodic')
close(fH)
% Anodic-leading
fH = figure();
imagesc(squeeze(mean(sinkAnodNthr(:,:,:)))')
axis square
colormap(flipud(gray()))
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\CSD_Ano_SinkMatrix_GrandMean_anodic')
close(fH)
% Difference
fH = figure();
imagesc(squeeze(mean(dif(:,:,:)))')
axis square
colormap(bluewhitered())
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\CSD_Ano_SinkMatrix_GrandMean_dif')
close(fH)

% Plot colorbars
% Gray
fH = figure();
colormap(flipud(gray()))
caxis([0 1])
colorbar()
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\CSD_Ano_SinkMatrix_colorbar_gray')
close(fH)
% Red-White-Blue
fH = figure();
caxis([-1 1])
colormap(bluewhitered())
colorbar()
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\CSD_Ano_SinkMatrix_colorbar_bwr')
close(fH)

% Calculate mean number of sinks per matrix
noSinksCath = mean(sum(sinkCathNthr,3),2);
noSinksAnod = mean(sum(sinkAnodNthr,3),2);
fprintf('Cathodic sink number: %4.3f +- %5.4f\n',mean(noSinksCath),std(noSinksCath))
fprintf('Anodic sink number: %4.3f +- %5.4f\n',mean(noSinksAnod),std(noSinksAnod))

% Figure S6A
fH = figure();
hold on
scatter(ones(1,8),noSinksCath,'filled','MarkerFaceColor',[0 0 120/255],'MarkerFaceAlpha',0.4)
scatter(ones(1,8)*2,noSinksAnod,'filled','MarkerFaceColor',[120/255 0 0],'MarkerFaceAlpha',0.4)
boxplot([noSinksCath,noSinksAnod])
set(gca,'YLim',[0 1.8])
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\CSD_Ano_SinkMatrix_noofsinks')
close(fH)

% Statistics
% [p,h,stats] = ranksum(noSinksCath,noSinksAnod)

%% Electric stimulation - Varying depth - CSD - Activity matrix differences
%
% Similiarity between Cathodic and Anodic activity matrices
%
% Euclidean distance between cathodic matrices calculated using normalized
% Frobenius norms
%
% Figure S6B
%
clear cathnorm anodnorm difnorm
cathnorm = NaN(8,8);
anodnorm = NaN(8,8);
for id1 = 8:-1:1
    for id2 = id1-1:-1:1
        cathnorm(id1,id2) = norm(double(squeeze(sinkCathNthr(id1,:,:)))-double(squeeze(sinkCathNthr(id2,:,:))),'fro')/norm(double(squeeze(sinkCathNthr(id1,:,:))),'fro');
        anodnorm(id1,id2) = norm(double(squeeze(sinkAnodNthr(id1,:,:)))-double(squeeze(sinkAnodNthr(id2,:,:))),'fro')/norm(double(squeeze(sinkAnodNthr(id1,:,:))),'fro');
    end
    difnorm(id1) = norm(double(squeeze(sinkCathNthr(id1,:,:)))-double(squeeze(sinkAnodNthr(id1,:,:))),'fro')/norm(double(squeeze(sinkCathNthr(id1,:,:))),'fro');
end

% Cathodic-leading
fH = figure();
imagesc(cathnorm)
colormap(flipud(gray()))
axis square
caxis([0 1.5])
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\CSD_Ano_SinkMatrix_FrobNorm_cathodic')
close(fH)

% Anodic-leading
fH = figure();
imagesc(anodnorm)
colormap(flipud(gray()))
axis square
caxis([0 1.5])
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\CSD_Ano_SinkMatrix_FrobNorm_anodic')
close(fH)

% gray colorbar
fH = figure();
colormap(flipud(gray()))
caxis([0 1.5])
colorbar()
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\CSD_Ano_SinkMatrix_FrobNomr_colorbar')
close(fH)

% Variance across experiments vs variance within experiment (= across polarity)

% Remove NaN entries in the activity matrices
cathnorm = cathnorm(:);
cathnorm(isnan(cathnorm))=[];

anodnorm = anodnorm(:);
anodnorm(isnan(anodnorm))=[];

% Figure 6A
fH = figure();
hold on
scatter(ones(1,28),cathnorm,'filled','MarkerFaceColor',[0 0 120/255],'MarkerFaceAlpha',0.2)
scatter(ones(1,28)*2,anodnorm,'filled','MarkerFaceColor',[120/255 0 0],'MarkerFaceAlpha',0.2)
scatter(ones(1,8)*3,difnorm,'filled','MarkerFaceColor',[0 0 0/255],'MarkerFaceAlpha',0.2)
boxplot([cathnorm;anodnorm;difnorm'],[ones(28,1);ones(28,1)*2;ones(8,1)*3])
set(gca,'YLim',[0 1.8])
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\CSD_Ano_Sink_FrobNorm_Difference_boxplot')
close(fH)

% Statistics
[p,h,stats] = kruskalwallis([cathnorm;anodnorm;difnorm'],[ones(28,1);ones(28,1)*2;ones(8,1)*3]) %#ok check results manually
% p > 0.05 -> the difference between activity matrices from different
% experiments is the same as the difference between cathodic/anodic 
% activity matrices of the same experiment

% Compare 1st quadrant with 4th quadrant
% Calculate nominal difference
clear cathnorm1 cathnorm2 anodnorm1 anodnorm2 difnorm1 difnorm2
cathnorm1 = NaN(8,8);
cathnorm2 = NaN(8,8);
anodnorm1 = NaN(8,8);
anodnorm2 = NaN(8,8);
for id1 = 8:-1:1
    for id2 = id1-1:-1:1
        cathnorm1(id1,id2) = norm(double(squeeze(sinkCathNthr(id1,1:8,1:8)))-double(squeeze(sinkCathNthr(id2,1:8,1:8))),'fro')/norm(double(squeeze(sinkCathNthr(id1,1:8,1:8))),'fro');
        cathnorm2(id1,id2) = norm(double(squeeze(sinkCathNthr(id1,9:16,9:16)))-double(squeeze(sinkCathNthr(id2,9:16,9:16))),'fro')/norm(double(squeeze(sinkCathNthr(id1,9:16,9:16))),'fro');
        anodnorm1(id1,id2) = norm(double(squeeze(sinkAnodNthr(id1,1:8,1:8)))-double(squeeze(sinkAnodNthr(id2,1:8,1:8))),'fro')/norm(double(squeeze(sinkAnodNthr(id1,1:8,1:8))),'fro');
        anodnorm2(id1,id2) = norm(double(squeeze(sinkAnodNthr(id1,9:16,9:16)))-double(squeeze(sinkAnodNthr(id2,9:16,9:16))),'fro')/norm(double(squeeze(sinkAnodNthr(id1,9:16,9:16))),'fro');
    end
    difnorm1(id1) = norm(double(squeeze(sinkCathNthr(id1,1:8,1:8)))-double(squeeze(sinkAnodNthr(id1,1:8,1:8))),'fro')/norm(double(squeeze(sinkCathNthr(id1,1:8,1:8))),'fro');
    difnorm2(id1) = norm(double(squeeze(sinkCathNthr(id1,9:16,9:16)))-double(squeeze(sinkAnodNthr(id1,9:16,9:16))),'fro')/norm(double(squeeze(sinkCathNthr(id1,9:16,9:16))),'fro');
end

% Figure 6B
fH = figure();
hold on
scatter(cathnorm1(:),cathnorm2(:),'filled','MarkerFaceColor',[0 0 120/255],'MarkerEdgeColor','none')
scatter(anodnorm1(:),anodnorm2(:),'filled','MarkerFaceColor',[120/255 0 0],'MarkerEdgeColor','none')
scatter(difnorm1,difnorm2,'filled','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none')
plot([0 2],[0 2],'--k')
axis square
set(gca,'XLim',[0 2],'YLim',[0 2])
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\CSD_Ano_Sink_FrobNorm_1stvs2ndhalf_scatter')
close(fH)

% Figure 6C - Cathodic-leading
fH = figure();
hold on
boxplot([cathnorm1(:),cathnorm2(:)],'Colors','k')
scatter(ones(1,64),cathnorm1(:),'filled','MarkerFaceColor',[0 0 120/255],'MarkerFaceAlpha',0.2)
scatter(ones(1,64)*2,cathnorm2(:),'filled','MarkerFaceColor',[0 0 120/255],'MarkerFaceAlpha',0.2)
set(gca,'YLim',[0 2])
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\CSD_Ano_Sink_FrobNorm_1stvs2ndhalf_boxCath')
close(fH)

% Figure 6C - Anodic-leading
fH = figure();
hold on
boxplot([anodnorm1(:),anodnorm2(:)],'Colors','k')
scatter(ones(1,64),anodnorm1(:),'filled','MarkerFaceColor',[120/255 0 0],'MarkerFaceAlpha',0.2)
scatter(ones(1,64)*2,anodnorm2(:),'filled','MarkerFaceColor',[120/255 0 0],'MarkerFaceAlpha',0.2)
set(gca,'YLim',[0 2])
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\CSD_Ano_Sink_FrobNorm_1stvs2ndhalf_boxAnod')
close(fH)

% Figure 6C - Difference
fH = figure();
hold on
boxplot([difnorm1(:),difnorm2(:)],'Colors','k')
scatter(ones(1,8),difnorm1(:),'filled','MarkerFaceColor',[0 0 0/255],'MarkerFaceAlpha',0.2)
scatter(ones(1,8)*2,difnorm2(:),'filled','MarkerFaceColor',[0 0 0/255],'MarkerFaceAlpha',0.2)
set(gca,'YLim',[0 2])
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\CSD_Ano_Sink_FrobNorm_1stvs2ndhalf_boxDif')
close(fH)

% Statistics
ranksum(cathnorm1(:),cathnorm2(:))
ranksum(anodnorm1(:),anodnorm2(:))
ranksum(difnorm1,difnorm2)


%% Electric stimulation - Varying depth - MUA - Load data

% Load multi-unit activity data for stimulation with varying depth
[~,~,MUASpikesA,MUASpikesC] = Load_MUA_Electric_VarDepth(cfg);

% Save MUA data to file for convenience
% save('D:\mbv\temp\Manuscript_anodic\anocath_layering_MUAs','MUAA','MUAC','MUASpikesA','MUASpikesC','-v7.3')
% load('D:\mbv\temp\Manuscript_anodic\anocath_layering_MUAs','MUAA','MUAC','MUASpikesA','MUASpikesC')

 

%% Electric stimulation - Varying depth - MUA - Figures

%
% Multi-unit activity images
%

% Spike raster plots - Overlay
% Stimulation on electrode 9
% Figure 7A, 7B
for el = 9
    
    % Cathodic-leading
    fH=figure('Renderer','Painters');
    hold on
    for e = 1:8
        for chan = 1:16
            for sw = 1:30
                if chan ~= el
                    plotdatax = find(MUASpikesC(e,el,chan,sw,990:1320));
                    plotdatay = -1*(sw+30*chan)*ones(length(find(MUASpikesC(e,el,chan,sw,990:1320))),1);
                    scatter(gca,plotdatax,plotdatay,'.','MarkerEdgeColor','k','MarkerFaceColor','k')
                else
                    continue
                end
            end
            if chan==el
                patch([0 462 462 0],[(-chan*30)-30 (-chan*30)-30 -chan*30 -chan*30],[1 1 0],'FaceAlpha',0.2,'EdgeColor','none')
            end
            plot([0 462],[(-chan*30)-30 (-chan*30)-30],':k')
        end
    end
    patch([110 176 176 110],[-510 -510 0 0],[0 0 0],'FaceAlpha',0.2,'EdgeColor','none')
    set(gca,'XTick',0:22:880,'XTickLabel',-5:100,'XLim',[0 330],'YLim',[-17*30 0])
    print(fH,'-dsvg','-r1200',['D:\mbv\temp\Manuscript_anodic\plots\Layer_Spk_cathodic_El' num2str(el)])
    close(fH)
    
    % Anodic-leading
    fH=figure('Renderer','Painters');
    hold on
    for e = 1:8
        for chan = 1:16
            for sw = 1:30
                if chan ~= el
                    plotdatax = find(MUASpikesA(e,el,chan,sw,990:1320));
                    plotdatay = -1*(sw+30*chan)*ones(length(find(MUASpikesA(e,el,chan,sw,990:1320))),1);
                    scatter(gca,plotdatax,plotdatay,'.','MarkerEdgeColor','k','MarkerFaceColor','k')
                else
                    continue
                end
            end
            if chan == el
                patch([0 462 462 0],[(-chan*30)-30 (-chan*30)-30 -chan*30 -chan*30],[1 1 0],'FaceAlpha',0.2,'EdgeColor','none')
            end
            plot([0 462],[(-chan*30)-30 (-chan*30)-30],':k')
        end
    end
    patch([110 176 176 110],[-510 -510 0 0],[0 0 0],'FaceAlpha',0.2,'EdgeColor','none')
    set(gca,'XTick',0:22:880,'XTickLabel',-5:100,'XLim',[0 330],'YLim',[-17*30 0])
    print(fH,'-dsvg','-r1200',['D:\mbv\temp\Manuscript_anodic\plots\Layer_Spk_anodic_El' num2str(el)])
    close(fH)
    
end


% Single-experiments
% Stimulation on electrode 9
%
% Figure S7
for e = 1:8
    % Cathodic-leading
    fH=figure('Renderer','Painters');
    hold on
        for chan = 1:16
            for sw = 1:30
                if chan ~= 9
                    plotdatax = find(MUASpikesC(e,9,chan,sw,1078:2200));
                    plotdatay = -1*(sw+30*chan)*ones(length(find(MUASpikesC(e,9,chan,sw,1078:2200))),1);
                    scatter(gca,plotdatax,plotdatay,'.','MarkerEdgeColor','k','MarkerFaceColor','k')
                else
                    continue
                end
            end
            if chan==9
                patch([0 242 242 0],[(-chan*30)-30 (-chan*30)-30 -chan*30 -chan*30],[1 1 0],'FaceAlpha',0.2,'EdgeColor','none')
            end
            plot([0 242],[(-chan*30)-30 (-chan*30)-30],':k')
        end
    patch([22 88 88 22],[-510 -510 0 0],[0 0 0],'FaceAlpha',0.2,'EdgeColor','none')
    set(gca,'XTick',0:22:660,'XTickLabel',-1:100,'XLim',[0 242],'YLim',[-17*30 0])
    print(fH,'-dsvg','-r1200',['D:\mbv\temp\Manuscript_anodic\plots\Layer_Spk_cathodic_El9_' cfg.Anoexps{e}])
    close(fH)
    
    % Anodic-leading
    fH=figure('Renderer','Painters');
    hold on
        for chan = 1:16
            for sw = 1:30
                if chan ~= 9
                    plotdatax = find(MUASpikesA(e,9,chan,sw,1078:1540));
                    plotdatay = -1*(sw+30*chan)*ones(length(find(MUASpikesA(e,9,chan,sw,1078:1540))),1);
                    scatter(gca,plotdatax,plotdatay,'.','MarkerEdgeColor','k','MarkerFaceColor','k')
                else
                    continue
                end
            end
            if chan == 9
                patch([0 242 242 0],[(-chan*30)-30 (-chan*30)-30 -chan*30 -chan*30],[1 1 0],'FaceAlpha',0.2,'EdgeColor','none')
            end
            plot([0 242],[(-chan*30)-30 (-chan*30)-30],':k')
        end
    patch([22 88 88 22],[-510 -510 0 0],[0 0 0],'FaceAlpha',0.2,'EdgeColor','none')
    set(gca,'XTick',0:22:660,'XTickLabel',-1:100,'XLim',[0 242],'YLim',[-17*30 0])
    print(fH,'-dsvg','-r1200',['D:\mbv\temp\Manuscript_anodic\plots\Layer_Spk_anodic_El9_' cfg.Anoexps{e}])
    close(fH)
    
end

% Calculate PSTH
clear pa pc
for e = 1:8
    allTSC = [];
    allTSA = [];
    for el = 9
        for ch=1:16
            for sw=1:30
                allTSC = [allTSC; find(squeeze(MUASpikesC(e,el,ch,sw,990:1320)))]; %#ok can't preallocate...
                allTSA = [allTSA; find(squeeze(MUASpikesA(e,el,ch,sw,990:1320)))]; %#ok can't preallocate...
            end
        end
    end
    allTSC = allTSC/22000*1000;
    allTSA = allTSA/22000*1000;
    
    [pc(e,:),edges] = hist(allTSC,0:0.5:15); %#ok can't preallocate...
    [pa(e,:),~] = hist(allTSA,0:0.5:15); %#ok can't preallocate...
end

% Figure 7C
fH = figure();
hold on
plot(mean(pc),'k')
plot(mean(pa),'r')
plot(mean(pc)+(std(pc)/sqrt(8)),'k')
plot(mean(pa)+(std(pa)/sqrt(8)),'r')
plot(mean(pc)-(std(pc)/sqrt(8)),'k')
plot(mean(pa)-(std(pa)/sqrt(8)),'r')
patch([21 33 33 21],[0 0 40 40],[0 0 0],'FaceAlpha',0.2,'EdgeColor','none')
set(gca,'YLim',[0 40],'XLim',[1 61],'XTick',1:4:65,'XTickLabel',-5:20)
print(fH,'-dsvg','-r1200',['D:\mbv\temp\Manuscript_anodic\plots\Layer_Spk_comparison_El' num2str(9)])
close(fH)


%% Electric stimulation - Varying depth - MUA - Analysis

% Sum of spikes in the 7 ms following the blanking interval, over all
% trials and all 15 recording electrodes
ampCtotal = zeros(8,16);
ampAtotal = zeros(8,16);
for e = 1:8
    for el = 1:16
        ampCtotal(e,el) = sum(sum(sum(sum(MUASpikesC(e,el,:,:,1166:1320),5),4),3));
        ampAtotal(e,el) = sum(sum(sum(sum(MUASpikesA(e,el,:,:,1166:1320),5),4),3));
    end
end

% Figure 8A
fH=figure();
hold on
errorbar(mean(ampCtotal),std(ampCtotal)/sqrt(8),'o-','Color',[0 0 120/255])
errorbar(mean(ampAtotal),std(ampAtotal)/sqrt(8),'o--','Color',[120/255 0 0])
set(gca,'YLim',[0 450])
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\Layer_Spk_totalAmp_layer')
close(fH)

% Figure 8B
fH=figure();
hold on
plot(ones(1,8),mean(ampCtotal,2),'o','MarkerFaceColor',[0 0 120/255],'MarkerEdgeColor','none')
plot(ones(1,8)*2,mean(ampAtotal,2),'o','MarkerFaceColor',[120/255 0 0],'MarkerEdgeColor','none')
boxplot([mean(ampCtotal,2),mean(ampAtotal,2)])
set(gca,'YLim',[0 450])
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\Layer_Spk_totalAmp_mean')
close(fH)

% Statistics for boxplot
% [p,h,stats] = signrank(mean(ampCtotal,2),mean(ampAtotal,2))
fprintf('MUA amplitude - cathodic: %4.2f +- %4.2f spikes\n',mean(mean(ampCtotal,2)),std(mean(ampCtotal,2)))
fprintf('MUA amplitude -   anodic: %4.2f +- %4.2f spikes\n',mean(mean(ampAtotal,2)),std(mean(ampAtotal,2)))

% repeated measures ANOVA
anim = categorical(repmat((1:8)',16,1));
stimel = categorical(reshape(repmat(1:16,8,1),[],1));
t = table(anim,stimel,ampCtotal(:),ampAtotal(:),'VariableNames',{'animal','StimEl','dC','dA'});
within = table(categorical([1 2]',[1 2],{'cathodic','anodic'}),'VariableNames',{'polarity'});
rm = fitrm(t,'dC,dA ~ StimEl','WithinDesign',within);
ranovatbl = ranova(rm,'WithinModel','polarity');
multcompare(rm,'polarity');
c = multcompare(rm,'polarity','By','StimEl');

% Calculate spike index as ratio of number of spikes in first ms over
% total number of spikes
ampC = zeros(8,16);
ampA = zeros(8,16);
for e = 1:8
    for el = 1:16
        ampC(e,el) = squeeze(sum(sum(sum(sum(MUASpikesC(e,el,:,:,1166:1188),5),4),3)))/ampCtotal(e,el);
        ampA(e,el) = squeeze(sum(sum(sum(sum(MUASpikesA(e,el,:,:,1166:1188),5),4),3)))/ampAtotal(e,el);
    end
end

% Figure 8C
fH=figure();
hold on
errorbar(mean(ampC),std(ampC)/sqrt(8),'Color',[0 0 120/255])
errorbar(mean(ampA),std(ampA)/sqrt(8),'Color',[120/255 0 0])
set(gca,'YLim',[0 1])
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\Layer_Spk_1ms')
close(fH)

% Figure 8D
fH=figure();
hold on
plot(ones(1,8),mean(ampC,2),'ok')
plot(ones(1,8)*2,mean(ampA,2),'or')
boxplot([mean(ampC,2),mean(ampA,2)])
set(gca,'XLim',[0 3],'YLim',[0 1])
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\Layer_Spk_1ms_boxplot')
close(fH)
% Statistics for boxplot
% [p,h,stats] = signrank(mean(ampC,2),mean(ampA,2))
fprintf('SI - cathodic: %4.2f +- %4.2f spikes\n',mean(mean(ampC,2)),std(mean(ampC,2)))
fprintf('SI -   anodic: %4.2f +- %4.2f spikes\n',mean(mean(ampA,2)),std(mean(ampA,2)))
