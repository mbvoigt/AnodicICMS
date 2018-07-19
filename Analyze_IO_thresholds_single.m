%% Analyze stimulation thresholds
%
% Subfunction for Anodic_Cathodic_main
%
% Determines stimulation threshold values based on CSD data
%
function [thrC,thrA] = Analyze_IO_thresholds_single(csdAnod,csdCath,strengths)
% Peak (sink) amplitude - baseline
ampBA = squeeze(max(csdAnod(:,:,:,1:1100),[],4));
ampBC = squeeze(max(csdCath(:,:,:,1:1100),[],4));

% Peak (sink) amplitude - response
ampA = squeeze(max(csdAnod(:,:,:,1100:2200),[],4));
ampC = squeeze(max(csdCath(:,:,:,1100:2200),[],4));

% Threshold determined as first response bigger than 2*baseline amplitude
thrC = zeros(8,16);
thrA = zeros(8,16);
for e = 1:8
    for ch = 1:16
        if any(ampC(e,:,ch)<ampBC(e,:,ch)*2)
            thrC(e,ch) = find(ampC(e,:,ch)<ampBC(e,:,ch)*2,1,'last'); 
        end
        if any(ampA(e,:,ch)<ampBA(e,:,ch)*2)
            thrA(e,ch) = find(ampA(e,:,ch)<ampBA(e,:,ch)*2,1,'last');
        end
    end
end

% If no response was found remove value from results
thrC(thrC==0)=NaN;
thrA(thrA==0)=NaN;

fprintf('mean %5.4f +- %5.4f\n',mean(strengths(min(thrC,[],2)),2),std(strengths(min(thrC,[],2)),[],2))
fprintf('mean %5.4f +- %5.4f\n',mean(strengths(min(thrA,[],2)),2),std(strengths(min(thrA,[],2)),[],2))
end