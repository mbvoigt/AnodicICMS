%% Analyze input/output function for electric stimulation
%
% Subfunction for Anodic_Cathodic_main
%
function [betaC,betaA,RsquaredC,RsquaredA,betaMC,betaMA,RsquaredMC,RsquaredMA] = Analyze_IO_curves(~,ampA,ampC)

% Strengths saved for convenience here:
strengths = [0.144000000000000 0.286000000000000 0.640888888888889 1.10972222222222 1.70572222222222 2.04505555555556 2.49694444444445 2.88944444444444 3.40527777777778 3.88100000000000 4.43355555555556 5.02438888888889 5.59955555555556 15.0966666666667 36.7231666666667];

% Initialize variables
betaC = zeros(8,4);
betaA = zeros(8,4);
RC = zeros(8,15);
RA = zeros(8,15);
SStotA = zeros(8,15);
SStotC = zeros(8,15);
RsquaredC = zeros(1,8);
RsquaredA = zeros(1,8);

% Non-Linear fit - Single experiments
for e = 1:8
    beta0 = [nanmean(nanmean(ampC(e,:,:),3),2)/2,nanmean(nanmean(ampC(e,:,:),3),2)*2,0.1,5];
    [betaC(e,:),RC(e,:)] = nlinfit(strengths,nanmean(nanmean(ampC(e,1:15,:),3),1),@(b,x)(b(1)+(b(2)./(1+exp(-b(3)*(x-b(4)))))),beta0);
    beta0 = [nanmean(nanmean(ampA(e,:,:),3),2)/2,nanmean(nanmean(ampA(e,:,:),3),2)*2,0.1,5];
    [betaA(e,:),RA(e,:)] = nlinfit(strengths,nanmean(nanmean(ampA(e,1:15,:),3),1),@(b,x)(b(1)+(b(2)./(1+exp(-b(3)*(x-b(4)))))),beta0);
    % Goodness of fit
    SStotC(e,:) = nanmean(nanmean(ampC(e,1:15,:),3),1) - nanmean(nanmean(ampC(e,1:15,:),3),2);
    SStotA(e,:) = nanmean(nanmean(ampA(e,1:15,:),3),1) - nanmean(nanmean(ampA(e,1:15,:),3),2);
    RsquaredC(e) = 1-(sum(RC(e,:).^2)./sum(SStotC(e,:).^2));
    RsquaredA(e) = 1-(sum(RA(e,:).^2)./sum(SStotA(e,:).^2));
end
clear SStotC SStotA RC RA beta0 e

% Non-linear fit - Mean of experiments
% Cathodic
% beta0 = [nanmean(nanmean(nanmean(ampC(:,:,:),3),1),2)/2,nanmean(nanmean(nanmean(ampC(:,:,:),3),1),2)*2,0.1,50];
beta0 = [nanmean(nanmean(ampC(:,1,:),3),1),nanmean(nanmean(ampC(:,15,:),3),1),1,5];
[betaMC,RMC] = nlinfit(strengths,nanmean(nanmean(ampC(:,1:15,:),3),1),@(b,x)(b(1)+(b(2)./(1+exp(-b(3)*(x-b(4)))))),beta0);
SStotMC = nanmean(nanmean(ampC(:,1:15,:),3),1) - nanmean(nanmean(nanmean(ampC(:,1:15,:),3),1),2);
RsquaredMC = 1-(sum(RMC.^2)./sum(SStotMC(:).^2));

% Anodic
% beta0 = [nanmean(nanmean(nanmean(ampA(:,:,:),3),1),2)/2,nanmean(nanmean(nanmean(ampA(:,:,:),3),1),2)*2,0.1,50];
beta0 = [nanmean(nanmean(nanmean(ampA(:,:,:),3),1),2)/2,nanmean(nanmean(nanmean(ampA(:,:,:),3),1),2)*2,0.1,5];
[betaMA,RMA] = nlinfit(strengths,nanmean(nanmean(ampA(:,1:15,:),3),1),@(b,x)(b(1)+(b(2)./(1+exp(-b(3)*(x-b(4)))))),beta0);
SStotMA = nanmean(nanmean(ampA(:,1:15,:),3),1) - nanmean(nanmean(nanmean(ampA(:,1:15,:),3),1),2);
RsquaredMA = 1-(sum(RMA.^2)./sum(SStotMA(:).^2));
clear RMC RMA SStotMA SStotMC beta0
