%% Accuracy

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
Acc_fixed = [];
Acc_online = [];
for s = 1:Nsub
    
    acc = load([fpath,'\Accuracy\',SubName{s}]);
    Acc_fixed(s,:) = acc.Acc_fixed*100;
    Acc_online(s,:) = acc.Acc_online*100;
end
%%

Acc_online = [0.4667 	0.2667	0.2667	0.4667	0.1333	0.4667;
1	0.7333	0.7333	0.8667	1	0.8667;
0.8667	0.8667	0.8667	0.9333	0.8		1;
0.7333	0.4667	0.6	0.7333	0.4	0.6667]*100;

Acc_fixed = [0.6667	0.4667	0.3333	0.2	0.2	1;
1	0.8	0.8667	0.7333	0.9333	1;
1	0.7333	0.8	0.7333	0.6667	1;
0.6	0.5333	0.5333	0.6	0.4667	0.8
]*100;
%% plot
figure; hold on;
colors = {[34, 139, 34]./255;[30, 144, 255]./255;[220, 20, 60]./255;[64, 224, 208]./255};
for s= 1:Nsub
a(s) = plot(Acc_online(s,:)','Color',colors{s},'marker','o','LineWidth',2,'MarkerFaceColor',colors{s});
end

a(Nsub+1) = plot(mean(Acc_online)','Color','black','marker','o','LineWidth',4,'MarkerFaceColor','black');

plot([1.5 1.5],[0 105],'Color',[0.3 0.3 0.3],'LineStyle','--')
plot([5.5 5.5],[0 105],'Color',[0.3 0.3 0.3],'LineStyle','--')

xlim([0.5 6.5])
ylim([0 105])
set(gca,'XTick',1:6,'XTickLabel',{'Pre','Main 1','Main 2','Main 3','Main 4','Post (adaptive)'},'FontSize',15)
ylabel('Accuracy (%)')
legend(a,{'Pilot01','Pilot02','Pilot03','Pilot04','Average'},'Location','eastoutside')

%%
Names = {'Pilot01','Pilot02','Pilot03','Pilot04'};
Acc_fixed_ = Acc_fixed;
Acc_fixed_(:,1) = Acc_online(:,1);
figure;
for s = 1:Nsub
    subplot(2,2,s);
    hold on;
    a1 = plot(Acc_online(s,:)','Color',colors{s},'marker','o','LineWidth',2,'MarkerFaceColor',colors{s});
    a2 = plot(Acc_fixed_(s,:)','Color',colors{s}-0.05,'marker','o','LineWidth',2,'MarkerFaceColor',colors{s},'LineStyle','--');
    plot([1.5 1.5],[0 105],'Color',[0.3 0.3 0.3],'LineStyle','--')
plot([5.5 5.5],[0 105],'Color',[0.3 0.3 0.3],'LineStyle','--')

xlim([0.5 6.5])
ylim([0 105])
set(gca,'XTick',1:6,'XTickLabel',{'Pre','Main 1','Main 2','Main 3','Main 4','Post'},'FontSize',15)
ylabel('Accuracy (%)')
legend([a1 a2],{'Adaptive','Fixed'},'Location','eastoutside')
title(Names{s})
end

