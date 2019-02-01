%% Analyze stimulation thresholds
%
% Subfunction for Anodic_Cathodic_main
%
% Used to plot CSD-IO figures of electric stimulation on electrodes 1,9,
% and 16
%
% Used to plot figures: 4A, S3A, S3B
%
function Plot_IO_figures(cfg,ampC,ampA,ampClog,ampAlog,strengths)

% Mean of experiments - Non-linear fit - Figure
% Figure 4A
fH = figure();
hold on
scatter(strengths,nanmean(nanmean(ampClog(:,1:15,:),3),1),'filled','MarkerEdgeColor','none','MarkerFaceColor',[0 0 120/255],'MarkerFaceAlpha',0.4)
errorbar(strengths,nanmean(nanmean(ampClog(:,1:15,:),3),1),nanstd(nanmean(ampClog(:,1:15,:),3),[],1)/sqrt(8),'-','Color',[0 0 120/255],'Marker','none')
scatter(strengths,nanmean(nanmean(ampAlog(:,1:15,:),3),1),'filled','MarkerEdgeColor','none','MarkerFaceColor',[120/255 0 0],'MarkerFaceAlpha',0.4)
errorbar(strengths,nanmean(nanmean(ampAlog(:,1:15,:),3),1),nanstd(nanmean(ampAlog(:,1:15,:),3),[],1)/sqrt(8),'-','Color',[120/255 0 0],'Marker','none')
set(gca,'YLim',[-10 20],'XScale','log')
print(fH,'-dsvg','-r1200',['D:\mbv\temp\Manuscript_anodic\plots\IO_CSD_NonLinFit_scaled_log_' cfg.saveFileString])
close(fH)


% Single experiments - Figure
% Figure S3A
% Cathodic
fH = figure();
hold on
for e = 1:8
    scatter(strengths,nanmean(ampC(e,1:15,:),3),'filled','MarkerEdgeColor','w','MarkerFaceColor',[0 0 120/255])
    plot(strengths,nanmean(ampC(e,1:15,:),3),'Color',[0 0 120/255])
end
set(gca,'YLim',[-15 25],'Xscale','log') 
% YLim upper limit: 
%       Electrode 1 & 16:  8000
%       Electrode 9:      35000
print(fH,'-dsvg','-r1200',['D:\mbv\temp\Manuscript_anodic\plots\IO_CSD_Cathodic_singleExperiments_log_' cfg.saveFileString])
close(fH)

% Anodic
% Figure S3B
fH = figure();
hold on
for e = 1:8
    scatter(strengths,nanmean(ampA(e,1:15,:),3),'filled','MarkerEdgeColor','w','MarkerFaceColor',[120/255 0 0])
    plot(strengths,nanmean(ampA(e,1:15,:),3),'Color',[120/255 0 0])
end
set(gca,'YLim',[-15 25],'Xscale','log')
% YLim upper limit: 
%       Electrode 1 & 16:  8000
%       Electrode 9:      35000
print(fH,'-dsvg','-r1200',['D:\mbv\temp\Manuscript_anodic\plots\IO_CSD_Anodic_singleExperiments_log_' cfg.saveFileString])
close(fH)
