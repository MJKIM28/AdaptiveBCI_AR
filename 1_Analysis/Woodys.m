%%
function [ERPaligned, LAG,ERPhistory,Laghistory,CORs] = Woodys(ERPs,limits,threshold,Fs,Template,ERPoutrange)
% ERPs: (time x channel x trial)
% limits: limits of lag (ms)
% threhold: stopping threshold
% Template: initial templte (time x channel)
% ERPoutrange: EEGs before and after 'ERPs'
%              ERPoutrange{1}: before ERPs (time x channel x trial) 
%              ERPoutrange{2}: after ERPs (time x channel x trial) 


limits = limits*Fs/1000;
[Ntime,Nch,Ntrial] = size(ERPs);

if nargin < 5
    Template = mean(ERPs,3,'omitnan');
end

if nargin < 6
    OutRangeFront = zeros(limits,Nch,Ntrial);
    OutRangeBack = zeros(limits,Nch,Ntrial);
else
    OutRangeFront = ERPoutrange{1}(end-limits+1:end,:,:);
    OutRangeBack = ERPoutrange{2}(1:limits,:,:);
end




% 3Hz lpf
[b,a] = butter(4,15./(Fs/2),'low');
% 15Hz lpf
[b2,a2] = butter(4,[1 15]./(Fs/2),'bandpass');

ERPhistory = [];
Laghistory = [];

CORs = [];
crit  = 1;
t = 1;
% figure;
while crit >= threshold
    %     Template = filtfilt(b2,a2,Template);
    Template_old = Template;
    COR = []; LAG = [];
    ERPaligned = [];
    for tr = 1:Ntrial
        %-- cross correlation (every channel)
        R_filt = [];
        for ch = 1:Nch
            [R,lag] = xcorr(Template(:,ch),ERPs(:,ch,tr),'normalized');
            %             R_filt(:,ch) = filtfilt(b,a,R);
            R_filt(:,ch) = R;

        end
        %-- mean of cross correlation over channels
        R_filtMean = mean(R_filt,2);

        %-- find peak lag
        %         [pks,locs] = findpeaks(R_filtMean);
        %         locsInLag = locs - (Ntime+1);
        %         [~,minId] = min(abs(locsInLag));
        %         Lag = locsInLag(minId);
        %         COR(tr) = pks(minId);
        %         LAG(tr) = Lag;

        %-- find max
        center = Ntime;
        R_InRange = R_filtMean(center-limits:center+limits);
        lagInRange = lag(center-limits:center+limits);
        [maxval,maxloc] = max(R_InRange);
        COR(tr) = maxval;
        Lag = lagInRange(maxloc);
        LAG(tr) = Lag;

        %-- realign ERPs
        if Lag < 0
            ERPaligned(:,:,tr) = cat(1,ERPs(abs(Lag)+1:end,:,tr), OutRangeBack(1:abs(Lag),:,tr));
        else
            ERPaligned(:,:,tr) = cat(1,OutRangeFront(end-abs(Lag)+1:end,:,tr), ERPs(1:end-abs(Lag),:,tr));
        end

    end
    ERPhistory(:,:,:,t) = ERPaligned;
    Laghistory(:,t) = LAG;

%     clf
%     subplot(121); plot(mean(ERPs,3));
%     subplot(122); plot(mean(ERPaligned,3));
%     sgtitle(num2str(t))
%     pause;

    CORmean = mean(COR);
    CORs(t) = CORmean;
    if t >1
        crit = CORs(t) - CORs(t-1);
    else
        crit = CORs(t);
    end
    t  = t  +1 ;
    Template = mean(ERPaligned,3,'omitnan');
end
end
