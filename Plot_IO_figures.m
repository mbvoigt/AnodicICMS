%% Analyze stimulation thresholds
%
% Subfunction for Anodic_Cathodic_main
%
% Used to plot CSD-IO figures of electric stimulation on electrodes 1,9,
% and 16
%
% Used to plot figures: 4A, S3A, S3B
%
function Plot_IO_figures(cfg,ampC,ampA,strengths,betaMC,betaMA)

% Mean of experiments - Non-linear fit - Figure
% Figure 4A
fH = figure();
hold on
scatter(strengths,nanmean(nanmean(ampC(:,1:15,:),3),1),'filled','MarkerEdgeColor','none','MarkerFaceColor',[0 0 120/255],'MarkerFaceAlpha',0.4)
errorbar(strengths,nanmean(nanmean(ampC(:,1:15,:),3),1),nanstd(nanmean(ampC(:,1:15,:),3),[],1)/sqrt(8),'o','Color',[0 0 120/255],'Marker','none')
plot(linspace(strengths(1),strengths(end),1000),betaMC(1)+(betaMC(2)./(1+exp(-betaMC(3)*(linspace(strengths(1),strengths(end),1000)-betaMC(4))))),'k-','MarkerSize',6,'MarkerEdgeColor','w','MarkerFaceColor','k')
scatter(strengths,nanmean(nanmean(ampA(:,1:15,:),3),1),'filled','MarkerEdgeColor','none','MarkerFaceColor',[120/255 0 0],'MarkerFaceAlpha',0.4)
errorbar(strengths,nanmean(nanmean(ampA(:,1:15,:),3),1),nanstd(nanmean(ampA(:,1:15,:),3),[],1)/sqrt(8),'o','Color',[120/255 0 0],'Marker','none')
plot(linspace(strengths(1),strengths(end),1000),betaMA(1)+(betaMA(2)./(1+exp(-betaMA(3)*(linspace(strengths(1),strengths(end),1000)-betaMA(4))))),'r-','MarkerSize',6,'MarkerEdgeColor','w','MarkerFaceColor','k')
set(gca,'YLim',[0 20000],'XScale','log')
print(fH,'-dsvg','-r1200',['D:\mbv\temp\Manuscript_anodic\plots\IO_CSD_NonLinFit_scaled_' cfg.saveFileString])
close(fH)


% Single experiments - Figure
% Figure S3A
% Cathodic
fH = figure();
hold on
for e = 1:8
    scatter(strengths,mean(ampC(e,1:15,:),3),'filled','MarkerEdgeColor','w','MarkerFaceColor',[0 0 120/255])
    plot(strengths,mean(ampC(e,1:15,:),3),'Color',[0 0 120/255])
end
set(gca,'YLim',[0 8000],'Xscale','log') 
% YLim upper limit: 
%       Electrode 1 & 16:  8000
%       Electrode 9:      35000
print(fH,'-dsvg','-r1200',['D:\mbv\temp\Manuscript_anodic\plots\IO_CSD_Cathodic_singleExperiments_' cfg.saveFileString])
close(fH)

% Anodic
% Figure S3B
fH = figure();
hold on
for e = 1:8
    scatter(strengths,mean(ampA(e,1:15,:),3),'filled','MarkerEdgeColor','w','MarkerFaceColor',[120/255 0 0])
    plot(strengths,mean(ampA(e,1:15,:),3),'Color',[120/255 0 0])
end
set(gca,'YLim',[0 8000],'Xscale','log')
% YLim upper limit: 
%       Electrode 1 & 16:  8000
%       Electrode 9:      35000
print(fH,'-dsvg','-r1200',['D:\mbv\temp\Manuscript_anodic\plots\IO_CSD_Anodic_singleExperiments_' cfg.saveFileString])
close(fH)
