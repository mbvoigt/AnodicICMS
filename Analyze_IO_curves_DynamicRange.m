%% Analyze Input-Output relationship for acoustical & electrical stimulation
%
% Subfunction for Anodic_Cathodic_main
%
% Plots figures: 4A, S5A, S5B
%
function [DRA,DRC,DRClick] = Analyze_IO_curves_DynamicRange(cfg,strengths,ampA,ampC,amp_click)

% current intensity in dB rel. to lowest stimulation current used (= 0.14 µA)
str_log = 10*log10(strengths/strengths(1));

% Amplitude normalization
for e = 1:8
    ampA(e,:,:) = ampA(e,:,:)./max(mean(ampA(e,:,:),3),[],2);
    ampC(e,:,:) = ampC(e,:,:)./max(mean(ampC(e,:,:),3),[],2);
    amp_click(e,:,:) = amp_click(e,:,:)./max(mean(amp_click(e,:,:),3),[],2);
end

% Normalized non-linear fits
%
% Single experiments
betaC_log = zeros(8,4);
betaA_log = zeros(8,4);
beta_click_log = zeros(8,4);
for e = 1:8
    beta0 = [nanmean(nanmean(ampC(e,:,:),3),2)/2,nanmean(nanmean(ampC(e,:,:),3),2)*2,0.1,50];
    betaC_log(e,:) = nlinfit(str_log,nanmean(nanmean(ampC(e,1:15,:),3),1),@(b,x)(b(1)+(b(2)./(1+exp(-b(3)*(x-b(4)))))),beta0);
    beta0 = [nanmean(nanmean(ampA(e,:,:),3),2)/2,nanmean(nanmean(ampA(e,:,:),3),2)*2,0.1,50];
    betaA_log(e,:) = nlinfit(str_log,nanmean(nanmean(ampA(e,1:15,:),3),1),@(b,x)(b(1)+(b(2)./(1+exp(-b(3)*(x-b(4)))))),beta0);
    beta0 = [nanmean(nanmean(amp_click(e,:,:),3),2)/2,nanmean(nanmean(amp_click(e,:,:),3),2)*2,0.1,50];
    beta_click_log(e,:) = nlinfit(1:17,nanmean(nanmean(amp_click(e,1:17,:),3),1),@(b,x)(b(1)+(b(2)./(1+exp(-b(3)*(x-b(4)))))),beta0);
end
clear beta0 e

% Population mean
%
% Cathodic-leading
beta0 = [nanmean(nanmean(nanmean(ampC(:,:,:),3),1),2)/2,nanmean(nanmean(nanmean(ampC(:,:,:),3),1),2)*2,0.1,50];
betaMC_log = nlinfit(str_log,nanmean(nanmean(ampC(:,1:15,:),3),1),@(b,x)(b(1)+(b(2)./(1+exp(-b(3)*(x-b(4)))))),beta0);
% Anodic-leading
beta0 = [nanmean(nanmean(nanmean(ampA(:,:,:),3),1),2)/2,nanmean(nanmean(nanmean(ampA(:,:,:),3),1),2)*2,0.1,50];
betaMA_log = nlinfit(str_log,nanmean(nanmean(ampA(:,1:15,:),3),1),@(b,x)(b(1)+(b(2)./(1+exp(-b(3)*(x-b(4)))))),beta0);
% Acoustic
beta0 = [nanmean(nanmean(nanmean(amp_click(:,:,:),3),1),2)/2,nanmean(nanmean(nanmean(amp_click(:,:,:),3),1),2)*2,0.1,50];
beta_clickM_log = nlinfit(1:17,nanmean(nanmean(amp_click(:,1:17,:),3),1),@(b,x)(b(1)+(b(2)./(1+exp(-b(3)*(x-b(4)))))),beta0);

% Align to a response amplitude equalling 0.5
syms x
eqn = betaMC_log(1)+((betaMC_log(2)-betaMC_log(1))./(1+exp(-betaMC_log(3)*(x-betaMC_log(4)))))==0.5;
offsetC = eval(solve(eqn,x));
eqn = betaMA_log(1)+((betaMA_log(2)-betaMA_log(1))./(1+exp(-betaMA_log(3)*(x-betaMA_log(4)))))==0.5;
offsetA = eval(solve(eqn,x));
eqn = beta_clickM_log(1)+((beta_clickM_log(2)-beta_clickM_log(1))./(1+exp(-beta_clickM_log(3)*(x-beta_clickM_log(4)))))==0.5;
offsetac = eval(solve(eqn,x));
clear x eqn

% Figure 4A, S5A, S5B
fH = figure();
hold on
scatter(str_log-offsetC,nanmean(nanmean(ampC(:,1:15,:),3),1),'filled','MarkerEdgeColor','none','MarkerFaceColor',[0 0 120/255],'MarkerFaceAlpha',0.4)
plot(linspace(str_log(1),str_log(end),1000)-offsetC,betaMC_log(1)+(betaMC_log(2)./(1+exp(-betaMC_log(3)*(linspace(str_log(1),str_log(end),1000)-betaMC_log(4))))),'k-','MarkerSize',6,'MarkerEdgeColor','w','MarkerFaceColor','k')
scatter(str_log-offsetA,nanmean(nanmean(ampA(:,1:15,:),3),1),'filled','MarkerEdgeColor','none','MarkerFaceColor',[120/255 0 0],'MarkerFaceAlpha',0.4)
plot(linspace(str_log(1),str_log(end),1000)-offsetA,betaMA_log(1)+(betaMA_log(2)./(1+exp(-betaMA_log(3)*(linspace(str_log(1),str_log(end),1000)-betaMA_log(4))))),'r-','MarkerSize',6,'MarkerEdgeColor','w','MarkerFaceColor','k')
scatter(((1:17)-offsetac)*5,nanmean(nanmean(amp_click(:,1:17,:),3),1),'filled','MarkerEdgeColor','none','MarkerFaceColor',[0 0 0],'MarkerFaceAlpha',0.4)
plot(((1:17)-offsetac)*5,beta_clickM_log(1)+(beta_clickM_log(2)./(1+exp(-beta_clickM_log(3)*((1:17)-beta_clickM_log(4))))),'k-','MarkerSize',6,'MarkerEdgeColor','w','MarkerFaceColor','k')
set(gca,'YLim',[0 1],'XLim',[-45 45])
print(fH,'-dsvg','-r1200',['D:\mbv\temp\Manuscript_anodic\plots\IO_CSD_AcousticvsElectric_' cfg.saveFileString])
close(fH)


%---------------------------
% Dynamic range calculation
%---------------------------

syms x
lowerC = zeros(8,1);
lowerA = zeros(8,1);
lowerClick = zeros(8,1);
upperC = zeros(8,1);
upperA = zeros(8,1);
upperClick = zeros(8,1);
for e = 1:8
% Lower limit -> 0.25 * max amplitude
eqn = betaC_log(e,1)+(betaC_log(e,2)./(1+exp(-betaC_log(e,3)*(x-betaC_log(e,4)))))==0.25;
lowerC(e) = eval(solve(eqn,x));
eqn = betaA_log(e,1)+(betaA_log(e,2)./(1+exp(-betaA_log(e,3)*(x-betaA_log(e,4)))))==0.25;
lowerA(e) = eval(solve(eqn,x));
eqn = beta_click_log(e,1)+(beta_click_log(e,2)./(1+exp(-beta_click_log(e,3)*(x-beta_click_log(e,4)))))==0.25;
lowerClick(e) = eval(solve(eqn,x));

% Upper limit -> 0.75 * max amplitude
eqn = betaC_log(e,1)+(betaC_log(e,2)./(1+exp(-betaC_log(e,3)*(x-betaC_log(e,4)))))==0.75;
upperC(e) = eval(solve(eqn,x));
eqn = betaA_log(e,1)+(betaA_log(e,2)./(1+exp(-betaA_log(e,3)*(x-betaA_log(e,4)))))==0.75;
upperA(e) = eval(solve(eqn,x));
eqn = beta_click_log(e,1)+(beta_click_log(e,2)./(1+exp(-beta_click_log(e,3)*(x-beta_click_log(e,4)))))==0.75;
upperClick(e) = eval(solve(eqn,x));
end
clear eqn x

% Remove complex values
for e = 1:8
    if imag(lowerC(e))~=0 || imag(lowerA(e))~=0 || imag(lowerClick(e))~=0
        lowerC(e) = NaN;
        lowerA(e) = NaN;
        lowerClick(e) = NaN;
        upperC(e) = NaN;
        upperA(e) = NaN;
        upperClick(e) = NaN;
    elseif imag(upperC(e))~=0 || imag(upperA(e))~=0 || imag(upperClick(e))~=0
        lowerC(e) = NaN;
        lowerA(e) = NaN;
        lowerClick(e) = NaN;
        upperC(e) = NaN;
        upperA(e) = NaN;
        upperClick(e) = NaN;
    end
end

% Dynamic range: Upper limit - Lower limit
DRC = upperC - lowerC;
DRA = upperA - lowerA;
DRClick = (upperClick - lowerClick)*5; % times five because of 5 dB steps

