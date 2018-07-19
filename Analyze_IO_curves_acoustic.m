%% Analyze Input-Output relationship for acoustical stimulation
%
% Subfunction for Anodic_Cathodic_main
%
% Plots figures: 4C, 4D
%
function [beta_click, beta_clickM, Rsquared_click] = Analyze_IO_curves_acoustic(cfg,amp_click,plotFigures)

if nargin < 3
    plotFigures = 0;
end

% Non-linear fit - Single experiments
beta_click = zeros(8,4); % parameters of non-linear fits
R = zeros(8,17); % residuals of non-linear fits 
SStot = zeros(8,17); % total sum of squares
Rsquared_click = zeros(8,1); % R^2
for e = 1:8
    fprintf('%s\n',cfg.SSAnoexps{e})
    beta0 = [mean(amp_click(e,1,:),3),mean(amp_click(e,17,:),3),0.5,10];
    [beta_click(e,:),R(e,:)] = nlinfit(1:17,mean(mean(amp_click(e,1:17,:),3),1),@(b,x)(b(1)+(b(2)./(1+exp(-b(3)*(x-b(4)))))),beta0);
    SStot(e,:) = mean(mean(amp_click(e,1:17,:),3),1) - mean(mean(amp_click(e,1:17,:),3),2);
    % Goodness of fit
    Rsquared_click(e) = 1-(sum(R(e,:).^2)./sum(SStot(e,:).^2));
end
clear beta0 e R

% Non-linear fit - Population mean
beta0 = [500,2000,0.5,8];
[beta_clickM,RM] = nlinfit(1:17,mean(mean(amp_click(:,1:17,:),3),1),@(b,x)(b(1)+(b(2)./(1+exp(-b(3)*(x-b(4)))))),beta0);
SStot = mean(mean(amp_click(:,1:17,:),3),1) - mean(mean(mean(amp_click(:,1:17,:),3),1),2);
% Goodness of fit
Rsquared_clickM = 1-(sum(RM.^2)./sum(SStot.^2)); %#ok value is actually not used
clear RM SStot beta0

%
% Figures
%
if plotFigures
    
% Single experiments
% Figure 4C 
fH = figure();
hold on
for e = 1:8
scatter(1:17,mean(mean(amp_click(e,1:17,:),3),1),'filled','MarkerEdgeColor','w','MarkerFaceColor','k')
plot(1:17,mean(amp_click(:,1:17,:),3),'k')
end
set(gca,'XTick',1:17,'XTickLabel',15:5:95)
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\IO_CSD_Acoustic_singleExperiments')
close(fH)

% Population mean
% Figure 4D
fH = figure();
hold on
plot(1:17,mean(mean(amp_click(:,1:17,:),3),1),'ok','MarkerEdgeColor','none','MarkerFaceColor','k')
plot(1:17,beta_clickM(1)+(beta_clickM(2)./(1+exp(-beta_clickM(3)*((1:17)-beta_clickM(4))))),'k-')
set(gca,'XTick',1:17,'XTickLabel',15:5:95)
print(fH,'-dsvg','-r1200','D:\mbv\temp\Manuscript_anodic\plots\IO_CSD_Acoustic_ExperimentMean')
close(fH)

clear fH
end
