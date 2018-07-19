%% Stimulation current calibration for electric stimulation
%
% Subfunction for Anodic_Cathodic_main
%
% Plots figures: S1A, S1B, S1C, S1D, S1E
%
function strengths = StimulationCurrentCalibration(StimCurrent_Cathodic,StimCurrent_Anodic)

% Current settings in AlphaOmega [µA]
strengths = [2,4,8,12,16,18,20,22,24,26,28,30,32,64,128];

% Remove experiments without current monitor
StimCurrent_Cathodic = StimCurrent_Cathodic([1,2,7],:,:,:);
StimCurrent_Anodic = StimCurrent_Anodic([1,2,7],:,:,:);

% Example trace
% Figure S1A left
fH = figure();
plot(squeeze(StimCurrent_Cathodic(1,13,10,47*44:53*44)),'Color',[0 0 120/255])
hold on
plot([3*44 3*44],[-12 12],'k:')
plot([3.2*44 3.2*44],[-12 12],'k:')
plot([3.4*44 3.4*44],[-12 12],'k:')
set(gca,'YLim',[-12 12],'XTick',0:44:440,'XTickLabel',-3:2)
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\Current_trace_cathodic')
close(fH)
% Figure S1A right
fH = figure();
plot(squeeze(StimCurrent_Anodic(3,13,10,47*44:53*44))','Color',[120/255 0 0])
hold on
plot([3*44 3*44],[-12 12],'k:')
plot([3.2*44 3.2*44],[-12 12],'k:')
plot([3.4*44 3.4*44],[-12 12],'k:')
set(gca,'YLim',[-12 12],'XTick',0:44:440,'XTickLabel',-3:2)
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\Current_trace_anodic')
close(fH)

% Current determination as peak amplitude of first phase:
strengthsC = abs(min(StimCurrent_Cathodic,[],4)); % absolute value of negative peak
strengthsA = max(StimCurrent_Anodic,[],4);  % peak value of positive phase

% Measured current vs amplitude setting
% Figure S1B left
% Cathodic
fH = figure();
hold on
for e = 1:3
    for sw = 1:30
        scatter(strengths,strengthsC(e,:,sw),'filled','MarkerFaceColor',[.8 .8 .8],'MarkerFaceAlpha',1,'SizeData',40)
    end
end
errorbar(strengths,mean(mean(strengthsC,3),1),std(mean(strengthsC,3),1),'Color',[0 0 120/255],'LineWidth',2)
set(gca,'YLim',[0 45])
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\Current_measure_cathodic')
close(fH)

% Figure S1B right
% Anodic
fH = figure();
hold on
for e = 1:3
    for sw = 1:30
        scatter(strengths,strengthsA(e,:,sw),'filled','MarkerFaceColor',[.8 .8 .8],'MarkerFaceAlpha',1,'SizeData',40)
    end
end
errorbar(strengths,mean(mean(strengthsA,3),1),std(mean(strengthsC,3),1),'Color',[120/255 0 0],'LineWidth',2)
set(gca,'YLim',[0 45])
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\Current_measure_anodic')
close(fH)

% Scatter plot: cathodic-vs-anodic
% Figure S1E
fH = figure();
scatter(strengthsC(:),strengthsA(:),'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.2,'SizeData',40)
hold on
plot([0 60],[0 60],':k')
set(gca,'XLim',[0 60],'YLim',[0 60])
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\Current_cathodic_vs_anodic')
close(fH)

% Difference index: Cathodic currents vs Anodic currents
for e = 1:3
    for str = 1:15
        for sw = 1:30
            strengthsD(e,str,sw) = (strengthsC(e,str,sw)-strengthsA(e,str,sw))./(strengthsC(e,str,sw)+strengthsA(e,str,sw));
        end
    end
end

% Mean over stimulus repetitions
strmean = mean(strengthsD,3);

% Difference index vs current setting
% Figure S1C left
fH = figure();
hold on
scatter(strengths,strmean(1,:),'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.2)
scatter(strengths,strmean(2,:),'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.2)
scatter(strengths,strmean(3,:),'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.2)
errorbar(strengths,mean(strmean),std(strmean),'k')
plot([0 140],[0 0],':k')
set(gca,'YLim',[-1 1])
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\Current_diff_currentsetting')
close(fH)

% Distribution
% Figure S1C middle
fH = figure();
scatter(zeros(1,45),strmean(:),'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.2)
hold on
plot([-0.5 0.5],[0 0],':k')
set(gca,'YLim',[-1 1])
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\Current_diff_distr')
close(fH)

% Histogram
% Figure S1C right
fH = figure();
histogram(strmean,20,'FaceColor','k','EdgeColor','none','FaceAlpha',1)
set(gca,'XLim',[-1 1])
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\Current_diff_hist')
close(fH)

% Variance over stimulus repetitions vs variance over experiments

% Cathodic-leading
fano_C_trial = squeeze(nanvar(strengthsC,[],3)./nanmean(strengthsC,3));
fano_C_exp = squeeze(nanvar(strengthsC,[],1)./nanmean(strengthsC,1));

% Figure S1D left
fH = figure();
bar(1,mean(mean(fano_C_trial)),'FaceColor',[0 0 120/255],'EdgeColor','none')
hold on
errorbar(1,mean(mean(fano_C_trial)),std(mean(fano_C_trial)),'k')
bar(2,mean(mean(fano_C_exp)),'FaceColor',[0 0 120/255],'EdgeColor','none')
errorbar(2,mean(mean(fano_C_exp)),std(mean(fano_C_exp)),'k')
set(gca,'YLim',[-0.1 0.4])
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\Current_fano_cathodic')
close(fH)

% Statistics
% [p,~,stats] = ranksum(mean(fano_C_trial),mean(fano_C_exp))

% Anodic-leading
fano_A_trial = squeeze(nanvar(strengthsA,[],3)./nanmean(strengthsA,3));
fano_A_exp = squeeze(nanvar(strengthsA,[],1)./nanmean(strengthsA,1));

% Figure S1D right
fH = figure();
bar(1,mean(mean(fano_A_trial)),'FaceColor','r','EdgeColor','none')
hold on
errorbar(1,mean(mean(fano_A_trial)),std(mean(fano_A_trial)),'k')
bar(2,mean(mean(fano_A_exp)),'FaceColor','r','EdgeColor','none')
errorbar(2,mean(mean(fano_A_exp)),std(mean(fano_A_exp)),'k')
set(gca,'YLim',[-0.1 0.4])
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\Current_fano_anodic')
close(fH)

% Statistics
% [p,~,stats] = ranksum(mean(fano_A_trial),mean(fano_A_exp))

% Actual strengths as means between cathodic- and anodic-leading
strengths = mean([mean(mean(strengthsC,3),1);mean(mean(strengthsA,3),1)],1);
