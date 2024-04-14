%% ERP comparison

clear all; close all

%% parameters
fpath = 'E:\1_연구\분석코드\2024_adaptive_BCI\Data\Pilot';

SubName = {'Subtest02','Subtest03','Subtest04','Subtest05'};
Nsub = length(SubName);

Fs = 500;
Nch = 31;
Tbase = 0.2;
Tepoc = 0.6;
Lepoc = (Tbase+Tepoc)*Fs;

Ntr_tr = 8;
Ntr_con  = 15;
Nsess_main = 4;

Nrepeat = 10;
chanloc = readlocs('D:\[1]EEGBCI\[2]Research\Code&Algorithm\BP_channelLocs.locs');

ChList = {'Fp1'    'Fpz'    'Fp2'    'F7'    'F3'    'Fz'    'F4'    'F8' ...
    'FT9'    'FC5'    'FC1'    'FC2'    'FC6'    'FT10'    'T7'    'C3'...
    'Cz'    'C4'    'T8'    'CP5'    'CP1'    'CP2'    'CP6'    'P7'...
    'P3'    'Pz'    'P4'    'P8'    'O1'    'Oz'    'O2'};
Times = -0.2:1/Fs:0.6-1/Fs;

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
    
    e_tr = load([fpath,'\Train\Epoch_15HZLPF\',SubName{s}]);
    e_te = load([fpath,'\Test\Epoch_15HZLPF\',SubName{s}]);
    chlist = setdiff(1:Nch,e_tr.Epoch.badch);

    %-- target
    EP_tar_tr(:,chlist,:,s) = e_tr.Epoch.tar;
    EP_tmp1 = [];
    for tr = 1:Ntr_con
        EP_tmp1 = cat(3,EP_tmp1,e_te.Epoch{tr}.tar);
    end
    EP_tmp2 = [];
    for tr = Ntr_con+1:Ntr_con*(Nsess_main+1)
        EP_tmp2 = cat(3,EP_tmp2,e_te.Epoch{tr}.tar);
    end
    EP_tmp3 = [];
    for tr = Ntr_con*(Nsess_main+1)+1:Ntr_con*(Nsess_main+2)
        EP_tmp3 = cat(3,EP_tmp3,e_te.Epoch{tr}.tar);
    end
    EP_tar_pre(:,chlist,:,s) = EP_tmp1;
    EP_tar_main(:,chlist,:,:,s) = reshape(EP_tmp2,[Lepoc,length(chlist),Nrepeat*Ntr_con,Nsess_main]);
    EP_tar_post(:,chlist,:,s) = EP_tmp3;   


    EP_nar_tr(:,chlist,:,s) = e_tr.Epoch.nar;
    EP_nar_tmp1 = [];
    for tr = 1:Ntr_con
        EP_nar_tmp1 = cat(3,EP_nar_tmp1,e_te.Epoch{tr}.nar(:,:,1:3*Nrepeat));
    end
    EP_nar_tmp2 = [];
    for tr = Ntr_con+1:Ntr_con*(Nsess_main+1)
        EP_nar_tmp2 = cat(3,EP_nar_tmp2,e_te.Epoch{tr}.nar);
    end
    EP_nar_tmp3 = [];
    for tr = Ntr_con*(Nsess_main+1)+1:Ntr_con*(Nsess_main+2)
        EP_nar_tmp3 = cat(3,EP_nar_tmp3,e_te.Epoch{tr}.nar);
    end
    EP_nar_pre(:,chlist,:,s) = EP_nar_tmp1;
    EP_nar_main(:,chlist,:,:,s) = reshape(EP_nar_tmp2,[Lepoc,length(chlist),Nrepeat*Ntr_con*3,Nsess_main]);
    EP_nar_post(:,chlist,:,s) = EP_nar_tmp3;  


end


ERPTarpre = squeeze(mean(EP_tar_pre,3,'omitnan'));
ERPTarmain = squeeze(mean(EP_tar_main,3,'omitnan'));
ERPTarpost = squeeze(mean(EP_tar_post,3,'omitnan'));

MeanTarpre = mean(ERPTarpre,3,'omitnan');
MeanTarmain = mean(ERPTarmain,4,'omitnan');
MeanTarpost = mean(ERPTarpost,3,'omitnan');


ERPNarpre = squeeze(mean(EP_nar_pre,3,'omitnan'));
ERPNarmain = squeeze(mean(EP_nar_main,3,'omitnan'));
ERPNarpost = squeeze(mean(EP_nar_post,3,'omitnan'));

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
color1 = {[0 0 0];[0.6 0.6 0.6]};
CondName = {'Pre','Post'};
%-- waveform
%-- individual
for s = 1:Nsub
    TarIn = cat(3,ERPDiffpre(:,:,s)',ERPDiffpost(:,:,s)');
topo_P300_eachcond_v2(ChList,TarIn,color1,CondName)
sgtitle(SubName{s})
end

%-- grand average
Tar = cat(3,MeanDiffpre',MeanDiffpost');
topo_P300_eachcond_v2(ChList,Tar,color1,CondName)

%-- amplitude


%%-- latency


%% main: 1 vs 2 vs 3 vs 4th session
color2 = {[0 0 0];[ 0, 100, 0]./255;[255, 165, 0]./255;[ 30, 144, 255]./255;[148, 0, 211]./255};
CondName2 = {'control','1st','2nd','3rd','4th'};
%-- individual
for s = 1:Nsub
topo_P300_eachcond_v2(ChList,permute(ERPDiffmain(:,:,:,s),[2,1,3]),color2,CondName2)
sgtitle(SubName{s})
end

%-- grand average
% topo_P300_eachcond_v2(ChList,permute(MeanDiffmain,[2,1,3]),color2,CondName2)
topo_P300_eachcond_v2(ChList,cat(3,MeanDiffpre',permute(MeanDiffmain,[2,1,3])),color2,CondName2)

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
CondName2 = {'Control','1st','2nd','3rd','4th'};

Wins = [0.1 0.2;0.2 0.3;0.3 0.4;0.4 0.5;0.5 0.6];
WinsName = {'100 ~ 200 ms';'200 ~ 300 ms';'300 ~ 400 ms';'400 ~ 500 ms';'500 ~ 600 ms';};
for tt = 1:size(Wins,1)
    for c = 1:Nsess_main+1

        subplot(5,5,(c-1)*length(Wins)+tt)
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

%% 4. latency (box plot)

%% Statistic
%-- Friedman
%% 1. amplitude

%% 2. latency




