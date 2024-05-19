%% ERP comparison

clear all; close all

%% parameters
fpath = 'E:\1_연구\분석코드\2024_adaptive_BCI\Data\Pilot';

SubNameList1 = {'Subtest02','Subtest03','Subtest04','Subtest05'};

SubNameList2 = {'Subtest12','Subtest13','Subtest14','Subtest15','Subtest17',...
    'Subtest18','Subtest19','Subtest20','Subtest21','Subtest22'};

SubName = SubNameList2;

Nsub = length(SubName);

Fs = 500;
Nch = 31;
Tbase = 0.2;
Tepoc = 0.6;
Lepoc = (Tbase+Tepoc)*Fs;

Ntr_tr = 15;
Ntr_con  = 15;
Nsess_main = 5;
Nsess = Nsess_main  + 2;
Nrepeat = 10;
chanloc = readlocs('D:\[1]EEGBCI\[2]Research\Code&Algorithm\BP_channelLocs.locs');

ChList = {'Fp1'    'Fpz'    'Fp2'    'F7'    'F3'    'Fz'    'F4'    'F8' ...
    'FT9'    'FC5'    'FC1'    'FC2'    'FC6'    'FT10'    'T7'    'C3'...
    'Cz'    'C4'    'T8'    'CP5'    'CP1'    'CP2'    'CP6'    'P7'...
    'P3'    'Pz'    'P4'    'P8'    'O1'    'Oz'    'O2'};
Times = -0.2:1/Fs:0.6-1/Fs;

CondName2 = {'Control','1st','2nd','3rd','4th','5th'};

load("TrialUsed.mat")
load("TrialUsed2.mat")
trial_use = [Trials(9:end); Trials2];
%% load (Epoch)
EP_tar_tr = NaN(Lepoc,Nch,Nrepeat*Ntr_tr,Nsub);
EP_tar_pre = NaN(Lepoc,Nch,Nrepeat*Ntr_con,Nsub);
EP_tar_main = NaN(Lepoc,Nch,Nrepeat*Ntr_con,Nsess_main,Nsub);
EP_tar_post = NaN(Lepoc,Nch,Nrepeat*Ntr_con,Nsub);
EP_nar_tr = NaN(Lepoc,Nch,Nrepeat*Ntr_tr*3,Nsub);
EP_nar_pre = NaN(Lepoc,Nch,Nrepeat*Ntr_con*3,Nsub);
EP_nar_main = NaN(Lepoc,Nch,Nrepeat*Ntr_con*3,Nsess_main,Nsub);
EP_nar_post = NaN(Lepoc,Nch,Nrepeat*Ntr_con*3,Nsub);
for s = 1:Nsub
     if s <5
        Ntr_con = 15;
    else
        Ntr_con = 14;
     end
    
    fprintf([SubName{s},'\n'])

    block_pre = trial_use{s}(1:find(trial_use{s} == Ntr_con));
    block_main = trial_use{s}(find(trial_use{s} == Ntr_con)+1:find(trial_use{s} == Ntr_con*(Nsess_main+1)));
    block_post = trial_use{s}(find(trial_use{s} == Ntr_con*(Nsess_main+1))+1:find(trial_use{s} == Ntr_con*(Nsess_main+2)));
    e_tr = load([fpath,'\Train\Epoch\',SubName{s}]);
    e_te = load([fpath,'\Test\Epoch\',SubName{s}]);
    chlist = setdiff(1:Nch,e_tr.Epoch.badch);

    %-- target
    EP_tar_tr(:,chlist,:,s) = e_tr.Epoch.tar;
    
    EP_tmp1 = [];
    for tr = 1:length(block_pre)
        EP_tmp1 = cat(3,EP_tmp1,e_te.Epoch{tr}.tar);
    end
    EP_tmp2 = [];
    for tr = length(block_pre)+1:length(block_pre)+length(block_main)
        EP_tmp2 = cat(3,EP_tmp2,e_te.Epoch{tr}.tar);
    end
    EP_tmp3 = [];
    for tr = length(block_pre)+length(block_main)+1:length(block_pre)+length(block_main)+length(block_post)
        EP_tmp3 = cat(3,EP_tmp3,e_te.Epoch{tr}.tar);
    end
    trials_pre = bsxfun(@plus,(block_pre-1)*Nrepeat,[1:Nrepeat]');
    trials_pre = trials_pre(:);
    EP_tar_pre(:,chlist,trials_pre,s) = EP_tmp1;

    trials_main = bsxfun(@plus,(block_main-Ntr_con-1)*Nrepeat,[1:Nrepeat]');
    trials_main = trials_main(:);
    EP_tmp2_ = NaN(Lepoc,Nch,Nrepeat*Ntr_con*Nsess_main);
    EP_tmp2_(:,chlist,trials_main) = EP_tmp2;
    EP_tar_main(:,:,1:Nrepeat*Ntr_con,:,s) = reshape(EP_tmp2_,[Lepoc,Nch,Nrepeat*Ntr_con,Nsess_main]);
    
    trials_post = bsxfun(@plus,(block_post-Ntr_con*(Nsess_main+1)-1)*Nrepeat,[1:Nrepeat]');
    trials_post = trials_post(:);
    EP_tar_post(:,chlist,trials_post,s) = EP_tmp3;


    EP_nar_tr(:,chlist,:,s) = e_tr.Epoch.nar;
    EP_nar_tmp1 = [];
    for tr = 1:length(block_pre)
        EP_nar_tmp1 = cat(3,EP_nar_tmp1,e_te.Epoch{tr}.nar(:,:,1:3*Nrepeat));
    end
    EP_nar_tmp2 = [];
    for tr = length(block_pre)+1:length(block_pre)+length(block_main)
        EP_nar_tmp2 = cat(3,EP_nar_tmp2,e_te.Epoch{tr}.nar);
    end
    EP_nar_tmp3 = [];
    for tr = length(block_pre)+length(block_main)+1:length(block_pre)+length(block_main)+length(block_post)
        EP_nar_tmp3 = cat(3,EP_nar_tmp3,e_te.Epoch{tr}.nar);
    end

    trials_preNt = bsxfun(@plus,(block_pre-1)*Nrepeat*3,[1:Nrepeat*3]');
    trials_preNt = trials_preNt(:);
    EP_nar_pre(:,chlist,trials_preNt,s) = EP_nar_tmp1;

    trials_mainNt = bsxfun(@plus,(block_main-Ntr_con-1)*Nrepeat*3,[1:Nrepeat*3]');
    trials_mainNt = trials_mainNt(:);
    EP_nar_tmp2_ = NaN(Lepoc,Nch,Nrepeat*3*Ntr_con*Nsess_main);
    EP_nar_tmp2_(:,chlist,trials_mainNt) = EP_nar_tmp2;
    EP_nar_main(:,:,1:Nrepeat*Ntr_con*3,:,s) = reshape(EP_nar_tmp2_,[Lepoc,Nch,Nrepeat*3*Ntr_con,Nsess_main]);
    
    trials_postNt = bsxfun(@plus,(block_post-Ntr_con*(Nsess_main+1)-1)*Nrepeat*3,[1:Nrepeat*3]');
    trials_postNt = trials_postNt(:);
    EP_nar_post(:,chlist,trials_postNt,s) = EP_nar_tmp3;



end

mkdir(['ERP'])
save(['ERP\ERP_tar'],'EP_tar_pre','EP_tar_post','EP_tar_main')
save(['ERP\ERP_nar'],'EP_nar_pre','EP_nar_post','EP_nar_main','-v7.3')

%%
load(['ERP\ERP_tar'])
load(['ERP\ERP_nar'])
sublist = 5:size(EP_tar_main,5);
ERPTarpre = squeeze(mean(EP_tar_pre(:,:,:,sublist),3,'omitnan'));
ERPTarmain = squeeze(mean(EP_tar_main(:,:,:,:,sublist),3,'omitnan'));
ERPTarpost = squeeze(mean(EP_tar_post(:,:,:,sublist),3,'omitnan'));

MeanTarpre = mean(ERPTarpre,3,'omitnan');
MeanTarmain = mean(ERPTarmain,4,'omitnan');
MeanTarpost = mean(ERPTarpost,3,'omitnan');


ERPNarpre = squeeze(mean(EP_nar_pre(:,:,:,sublist),3,'omitnan'));
ERPNarmain = squeeze(mean(EP_nar_main(:,:,:,:,sublist),3,'omitnan'));
ERPNarpost = squeeze(mean(EP_nar_post(:,:,:,sublist),3,'omitnan'));

MeanNarpre = mean(ERPNarpre,3,'omitnan');
MeanNarmain = mean(ERPNarmain,4,'omitnan');
MeanNarpost = mean(ERPNarpost,3,'omitnan');


ERPDiffpre = ERPTarpre - ERPNarpre;
ERPDiffmain = ERPTarmain - ERPNarmain;
ERPDiffpost = ERPTarpost - ERPNarpost;

MeanDiffpre = mean(ERPDiffpre,3,'omitnan');
MeanDiffmain = mean(ERPDiffmain,4,'omitnan');
MeanDiffpost = mean(ERPDiffpost,3,'omitnan');

%% ERP 
%% pre vs. post
Nsub = length(sublist);
color1 = {[0 0 0];[0.6 0.6 0.6]};
CondName = {'Pre','Post'};
%-- waveform
%-- individual
for s = 1:Nsub
    TarIn = cat(3,ERPDiffpre(:,:,s)',ERPDiffpost(:,:,s)');
topo_P300_eachcond_v2(ChList,TarIn,color1,CondName,[-8 8])
sgtitle(SubName{s})
end

%-- grand average
Tar = cat(3,MeanDiffpre',MeanDiffpost');
topo_P300_eachcond_v2(ChList,Tar,color1,CondName,[-3.5 3.5])

%% -- amplitude & latency
Win{1} = find(Times > 0.1 & Times < 0.3);
Win{2} = find(Times > 0.25 & Times < 0.45);
winsize = fix(0.025*Fs);

AmpTar = NaN(Nch,length(Win),2,Nsub,Nsess);
PeakLatTar = NaN(Nch,length(Win),2,Nsub,Nsess);
AmpNar = NaN(Nch,length(Win),2,Nsub,Nsess);
PeakLatNar = NaN(Nch,length(Win),2,Nsub,Nsess);
for s = 1:Nsub
    chlist = ~isnan(ERPDiffpre(1,:,s));
        X = ERPDiffpre(:,chlist,s)';
        [AmpTar(chlist,:,:,s,1),PeakLatTar(chlist,:,:,s,1)] = getPeakAmp(X,Win,winsize);
      
        for ss = 1:Nsess_main
         X = ERPDiffmain(:,chlist,ss,s)';
        [AmpTar(chlist,:,:,s,ss+1),PeakLatTar(chlist,:,:,s,ss+1)] = getPeakAmp(X,Win,winsize);
        end
      
         X = ERPDiffpost(:,chlist,s)';
        [AmpTar(chlist,:,:,s,end),PeakLatTar(chlist,:,:,s,end)] = getPeakAmp(X,Win,winsize);
      
end
save('ERP\PeakTar','AmpTar','PeakLatTar')

chdraw = find(ismember(ChList,{'Fz','FC1','FC2','Cz','P3','Pz','P4','O1','Oz','O2'}));
% P2
P2amp = squeeze(AmpTar(chdraw(1:4),1,2,:,:));
P2lat = squeeze(PeakLatTar(chdraw(1:4),1,2,:,:));
P2lat(isnan(P2lat)) = 1;
P2lat = Times(P2lat);
P2lat(P2lat<0) = NaN;

% P3
P3amp = squeeze(AmpTar(chdraw(5:7),2,2,:,:));
P3lat = squeeze(PeakLatTar(chdraw(5:7),2,2,:,:));
P3lat(isnan(P3lat)) = 1;
P3lat = Times(P3lat);
P3lat(P3lat<0) = NaN;

% N2
N2amp = squeeze(AmpTar(chdraw(8:10),1,1,:,:));
N2lat = squeeze(PeakLatTar(chdraw(8:10),1,1,:,:));
N2lat(isnan(N2lat)) = 1;
N2lat = Times(N2lat);
N2lat(N2lat<0) = NaN;

save(['ERP\value_ERPComp_FC1FC2'],'P2amp','P2lat','P3amp','P3lat',...
    'N2amp','N2lat')
%% main: 1 vs 2 vs 3 vs 4 vs 5th session
color2 = {[0 0 0];[0, 0, 255]./255;[255, 0, 0]./255;[0, 128, 0]./255;[255, 165, 0]./255;[128, 0, 128]./255};

%-- individual
for s = 1:Nsub
topo_P300_eachcond_v2(ChList,cat(3,ERPDiffpre(:,:,s)',permute(ERPDiffmain(:,:,:,s),[2,1,3])),color2,CondName2,[-5 5])
sgtitle(SubName{s})
end

%-- grand average
% topo_P300_eachcond_v2(ChList,permute(MeanDiffmain,[2,1,3]),color2,CondName2)
topo_P300_eachcond_v2(ChList,cat(3,MeanDiffpre',permute(MeanDiffmain,[2,1,3])),color2,CondName2,[-3.5 3.5])

%-- amplitude

%-- latency


%%
%% plot
%% 1. waveform
chdrawname = {'Fz','FC1','FC2','Cz','P3','Pz','P4','O1','Oz','O2'};
chdrawid = find(ismember(ChList,chdrawname));
topo_P300_partial(chdrawname,cat(3,MeanDiffpre(:,chdrawid,:)',permute(MeanDiffmain(:,chdrawid,:),[2,1,3])),color2,CondName2)

color1 = {[0 0 0];[0.6 0.6 0.6]};
CondName = {'Pre','Post'};
topo_P300_partial(chdrawname,cat(3,MeanDiffpre(:,chdrawid,:)',permute(MeanDiffpost(:,chdrawid,:),[2,1,3])),color1,CondName)

%% 2. Topography

%-- topoplot
figure;


Wins = [0.1 0.2;0.2 0.3;0.3 0.4;0.4 0.5;0.5 0.6];
WinsName = {'100 ~ 200 ms';'200 ~ 300 ms';'300 ~ 400 ms';'400 ~ 500 ms';'500 ~ 600 ms';};
for tt = 1:size(Wins,1)
    for c = 1:Nsess_main+1

        subplot(6,5,(c-1)*length(Wins)+tt)
        if c == 1
        X = squeeze(mean(MeanDiffpre(Times>Wins(tt,1) & Times<Wins(tt,2),:),'omitnan'));
        else
        X = squeeze(mean(MeanDiffmain(Times>Wins(tt,1) & Times<Wins(tt,2),:,c-1),'omitnan'));
        end
        topoplot(X,chanloc);
        caxis([-0.8 0.8])

        if c == 1
            title(WinsName{tt})
        end
        if tt == 1
            ylabel(CondName2{c})
        end

    end
end
set(gca,'FontSize',15)
%% 3. amplitude (box plot)
figure;
pb = [];
ch1 = 0; ch2 = 0; ch3 = 0;
color3 = [color2; [0 0 0]];
CondName3 = ['Pre',CondName2(2:end),'Post'];
for ch = 1:length(chdraw)
%     axs{ch} = subplot(4,3,chdrawID(ch));
if strcmp(chdrawname{ch},'Fz')
        axs{ch} = axes('Position',[0.4208    0.83    0.18    0.1430]);
    elseif strcmp(chdrawname{ch},'FC1')
        axs{ch} =axes('Position',[0.18 0.65  0.18   0.1430]);
    elseif strcmp(chdrawname{ch},'FC2')
       axs{ch} = axes('Position',[0.64 0.65  0.18   0.1430]);
    elseif strcmp(chdrawname{ch},'Cz')
       axs{ch} = axes('Position',[0.4108    0.55    0.18    0.1430]);
    elseif strcmp(chdrawname{ch},'P3')
       axs{ch} = axes('Position',[0.05    0.3291    0.18    0.1430]);
    elseif strcmp(chdrawname{ch},'Pz')
        axs{ch} =axes('Position',[0.4108    0.3291   0.18    0.1430]);
    elseif strcmp(chdrawname{ch},'P4')
       axs{ch} = axes('Position',[ 0.78    0.3291    0.18    0.1430]);
    elseif strcmp(chdrawname{ch},'O1')
       axs{ch} = axes('Position',[ 0.18    0.1100   0.18   0.1430]);
    elseif strcmp(chdrawname{ch},'Oz')
       axs{ch} = axes('Position',[ 0.4108    0.1100   0.18   0.1430]);
    elseif strcmp(chdrawname{ch},'O2')
       axs{ch} = axes('Position',[0.64    0.1100    0.18   0.1430]);
end

    hold on

    Chlabel = ChList{chdraw(ch)};

    if ch <= 4

        ch1 = ch1 + 1;
        Ys = squeeze(P2amp(ch1,:,:));

    elseif ch > 4 && ch <=7
        ch2 = ch2 + 1;
        Ys = squeeze(P3amp(ch2,:,:));
    elseif ch > 7 && ch <=10
        ch3 = ch3 + 1;
        Ys = squeeze(N2amp(ch3,:,:));
    end



    for c = 1:Nsess
        subind = ~isnan(Ys(:,c));
        errorbar(c,abs(mean(Ys(subind,c))),std(Ys(subind,c))/sqrt(sum(subind)),'linewidth',2,'color',color3{c},'Marker','square');
        %          boxchart(ones(Nsub,1)*c,Ys(:,c),'BoxFaceColor',colors{ConditionOrder(c)},'MarkerStyle','.','jitteroutliers','on','MarkerColor',colors{ConditionOrder(c)})
    end

    if ch == 1
        ylabel('Amplitude (\muV)')
        set(gca,'xtick',1:Nsess,'XTicklabel',CondName3,'XTickLabelRotation',45)
    else
        set(gca,'xtick',1:Nsess,'XTicklabel',[])
    end

    title(Chlabel)

    set(gca,'fontsize',12,'Box','off','TickDir','out')
    xlim([0.5 7.5])
    if ch <= 4
            ylim([0.3 2.5])
    elseif ch <= 7
         ylim([0 2.5])
    else
        ylim( [1 5])
    end
end
% 
% [~,~,~,p_ad] =fdr_bh(p)
% 
% for ch = 1:length(chdraw)
%     post{ch} = multcompare(stats{ch},"Display","off");
% end

%% 4. latency (box plot)
%% 3. amplitude (box plot)
figure;
pb = [];
ch1 = 0; ch2 = 0; ch3 = 0;
color3 = [color2; [0 0 0]];
CondName3 = ['Pre',CondName2(2:end),'Post'];
for ch = 1:length(chdraw)
%     axs{ch} = subplot(4,3,chdrawID(ch));
if strcmp(chdrawname{ch},'Fz')
        axs{ch} = axes('Position',[0.4208    0.83    0.18    0.1430]);
    elseif strcmp(chdrawname{ch},'FC1')
        axs{ch} =axes('Position',[0.18 0.65  0.18   0.1430]);
    elseif strcmp(chdrawname{ch},'FC2')
       axs{ch} = axes('Position',[0.64 0.65  0.18   0.1430]);
    elseif strcmp(chdrawname{ch},'Cz')
       axs{ch} = axes('Position',[0.4108    0.55    0.18    0.1430]);
    elseif strcmp(chdrawname{ch},'P3')
       axs{ch} = axes('Position',[0.05    0.3291    0.18    0.1430]);
    elseif strcmp(chdrawname{ch},'Pz')
        axs{ch} =axes('Position',[0.4108    0.3291   0.18    0.1430]);
    elseif strcmp(chdrawname{ch},'P4')
       axs{ch} = axes('Position',[ 0.78    0.3291    0.18    0.1430]);
    elseif strcmp(chdrawname{ch},'O1')
       axs{ch} = axes('Position',[ 0.18    0.1100   0.18   0.1430]);
    elseif strcmp(chdrawname{ch},'Oz')
       axs{ch} = axes('Position',[ 0.4108    0.1100   0.18   0.1430]);
    elseif strcmp(chdrawname{ch},'O2')
       axs{ch} = axes('Position',[0.64    0.1100    0.18   0.1430]);
end

    hold on

    Chlabel = ChList{chdraw(ch)};

    if ch <= 4

        ch1 = ch1 + 1;
        Ys = squeeze(P2lat(ch1,:,:));

    elseif ch > 4 && ch <=7
        ch2 = ch2 + 1;
        Ys = squeeze(P3lat(ch2,:,:));
    elseif ch > 7 && ch <=10
        ch3 = ch3 + 1;
        Ys = squeeze(N2lat(ch3,:,:));
    end

    for c = 1:Nsess
        subind = ~isnan(Ys(:,c));
        errorbar(c,abs(mean(Ys(subind,c))),std(Ys(subind,c))/sqrt(sum(subind)),'linewidth',2,'color',color3{c},'Marker','square');
        %          boxchart(ones(Nsub,1)*c,Ys(:,c),'BoxFaceColor',colors{ConditionOrder(c)},'MarkerStyle','.','jitteroutliers','on','MarkerColor',colors{ConditionOrder(c)})
    end

    if ch == 1
        ylabel('Latency (s)')
        set(gca,'xtick',1:Nsess,'XTicklabel',CondName3,'XTickLabelRotation',45)
    else
        set(gca,'xtick',1:Nsess,'XTicklabel',[])
    end

    title(Chlabel)

    set(gca,'fontsize',12,'Box','off','TickDir','out')
    xlim([0.5 7.5])
    if ch <= 4
            ylim([0.15 0.3])
    elseif ch <= 7
         ylim([0.25 0.45])
    else
        ylim( [0.17 0.27])
    end
end
set(gcf,'position',[263   227  1079  752])
%% Statistic
%-- Friedman
%% 1. amplitude

%% 2. latency




