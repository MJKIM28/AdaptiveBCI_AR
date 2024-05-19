



SubNames = {'Subtest12','Subtest13','Subtest14','Subtest15','Subtest17',...
    'Subtest18','Subtest19','Subtest20','Subtest21','Subtest22'};
Nsub = length(SubNames);

Nch = 31;
sublist = 5:Nsub;

Fs = 500;
Time2 = 0:1/Fs:0.6-1/Fs;

Windowlength = fix(0.2*Fs);
Windowlist = fix(1:0.05*Fs:0.6*Fs-Windowlength);
Ntime = length(Windowlist);

load('adaptiveSPHDCA\SS03_exp_adaptiveSPHDCA_LDA_l1_0.03_l2_0.03_thres1_1_thres2_1.mat')
chanloc = readlocs('D:\[1]EEGBCI\[2]Research\Code&Algorithm\BP_channelLocs.locs');
%%
DSP_init = NaN(Nch,Ntime,length(sublist));
DSP_init2 = NaN(Nch,Ntime,length(sublist));
DSP_init3 = NaN(Nch,Ntime,length(sublist));
DSP_tot = NaN(Nch,Nch,length(sublist));

ss = 1;
for s = 5:Nsub
    if s == 5
        p = 3;
    else
         p = 2;
    end
    chlist = setdiff(1:Nch,BadChs{s});
DSP_init(chlist,:,ss) = squeeze(DSPs{1,p,s}(:,1,:));
DSP_init2(chlist,:,ss) = squeeze(DSPs{1,p,s}(:,2,:));
DSP_init3(chlist,:,ss) = squeeze(DSPs{1,p,s}(:,3,:));
DSP_tot(chlist,chlist,ss) = squeeze(DSPtotal{s});

ss = ss +1;
end

DSPinv_init = NaN(Nch,Ntime,length(sublist));
ss = 1;
for s = 5:Nsub
    if s == 5
        p = 3;
    else
         p = 2;
    end
    chlist = setdiff(1:Nch,BadChs{s});
    
    W = [];
    for t = 1:Ntime
    W(:,:,t) = inv(DSPs{1,p,s}(:,:,t));
    end
DSPinv_init(chlist,:,ss) = squeeze(W(:,1,:));

ss = ss +1;
end


figure; 
topoplot(DSP_init(:,end,2),chanloc)
figure; 
topoplot(DSP_init2(:,end,2),chanloc)
figure; 
topoplot(DSP_tot(:,1,2),chanloc)

%%
figure;
for t = 1:Ntime
subplot(1,Ntime,t)
topoplot(mean(abs(DSP_init(:,t,:)),3,'omitnan'),chanloc,'headrad',0.62);
caxis([0 .5])
title([num2str(Time2(Windowlist(t))*1000),'ms ~ ',num2str(Time2(Windowlist(t))*1000+200),'ms'])
end
figure;
topoplot(mean(abs(DSP_tot(:,1,:)),3,'omitnan'),chanloc,'headrad',0.62);
caxis([0 .5])

figure;
for t = 1:Ntime
subplot(1,Ntime,t)
topoplot(mean(abs(DSPinv_init(:,t,:)),3,'omitnan'),chanloc)
% caxis([0 9])
title([num2str(Time2(Windowlist(t))*1000),'ms ~ ',num2str(Time2(Windowlist(t))*1000+200),'ms'])
colorbar
end