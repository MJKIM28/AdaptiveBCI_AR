%% Accuracy

clear all; close all
%%
filenames = {'Subtest07','Subtest08','Subtest09','Subtest10','Subtest12','Subtest13','Subtest14','Subtest15'}
% for ss = 1:length(filenames)
%     load(filenames{ss})
%     Acc_fixed = Acc_fix;
%     Acc_online = Acc_ad;
%     save(filenames{ss},'Acc_online','Acc_fixed');
% end

for ss = 6:length(filenames)
    load(filenames{ss})
    if length(Acc_online) < 7
    Acc_online = [Acc_fixed(1) Acc_online];
    end
    save(filenames{ss},'Acc_online','Acc_fixed');
end


%% parameters
fpath = 'E:\1_연구\분석코드\2024_adaptive_BCI\Data\Pilot';

SubName = {'Subtest02','Subtest03','Subtest04','Subtest05','Subtest07','Subtest08','Subtest09','Subtest10','Subtest12','Subtest13','Subtest14','Subtest15'};
Nsub = length(SubName);

Fs = 500;
Nch = 31;
Tbase = 0.2;
Tepoc = 0.6;
Lepoc = (Tbase+Tepoc)*Fs;

Ntr_tr = 8;
Ntr_con  = 15;
Nsess_main = 5;

Nrepeat = 10;
chanloc = readlocs('D:\[1]EEGBCI\[2]Research\Code&Algorithm\BP_channelLocs.locs');

ChList = {'Fp1'    'Fpz'    'Fp2'    'F7'    'F3'    'Fz'    'F4'    'F8' ...
    'FT9'    'FC5'    'FC1'    'FC2'    'FC6'    'FT10'    'T7'    'C3'...
    'Cz'    'C4'    'T8'    'CP5'    'CP1'    'CP2'    'CP6'    'P7'...
    'P3'    'Pz'    'P4'    'P8'    'O1'    'Oz'    'O2'};
Times = -0.2:1/Fs:0.6-1/Fs;


% Define the start and end colors
start_color = [0, 0, 255]; % Blue
end_color = [255, 0, 0]; % Red

% Number of samples
n_samples = 12;

% Create the gradient
gradient = zeros(Nsub, 3); % Initialize the matrix for RGB values
for i = 1:3
    gradient(:, i) = linspace(start_color(i), end_color(i), n_samples);
end

%% load (Epoch)
Acc_fixed = NaN(Nsub,7);
Acc_online = NaN(Nsub,7);
for s = 1:Nsub

    acc = load([fpath,'\Accuracy\',SubName{s}]);
    if s < 5
        Acc_fixed(s,[1:5,7]) = acc.Acc_fixed*100;
        Acc_online(s,[1:5,7]) = acc.Acc_online*100;
    else
        Acc_fixed(s,:) = acc.Acc_fixed*100;
        Acc_online(s,:) = acc.Acc_online*100;
    end
end
%%

% Acc_online = [0.4667 	0.2667	0.2667	0.4667	0.1333	0.4667;
% 1	0.7333	0.7333	0.8667	1	0.8667;
% 0.8667	0.8667	0.8667	0.9333	0.8		1;
% 0.7333	0.4667	0.6	0.7333	0.4	0.6667]*100;
% 
% Acc_fixed = [0.6667	0.4667	0.3333	0.2	0.2	1;
% 1	0.8	0.8667	0.7333	0.9333	1;
% 1	0.7333	0.8	0.7333	0.6667	1;
% 0.6	0.5333	0.5333	0.6	0.4667	0.8
% ]*100;
Acc_fixed(3,:) = [ 1	0.7333	0.8	0.7333	0.6667 NaN	1]*100
Acc_online(4,:) = [0.7333	0.4667	0.6	0.7333	0.4	0.6667]*100
%% plot
Names = {'Pilot01','Pilot02','Pilot03','Pilot04','Pilot06','Pilot07','Pilot08','Pilot09','Pilot11','Pilot12','Pilot13','Pilot14'};

figure; hold on;
colors = gradient./255;
for s= 1:Nsub
a(s) = plot(Acc_online(s,:)','Color',colors(s,:),'marker','o','LineWidth',2,'MarkerFaceColor',colors(s,:));
end

a(Nsub+1) = plot(mean(Acc_online,'omitnan')','Color','black','marker','o','LineWidth',4,'MarkerFaceColor','black');

plot([1.5 1.5],[0 105],'Color',[0.3 0.3 0.3],'LineStyle','--')
plot([6.5 6.5],[0 105],'Color',[0.3 0.3 0.3],'LineStyle','--')

xlim([0.5 7.5])
ylim([0 105])
set(gca,'XTick',1:7,'XTickLabel',{'Pre','Main 1','Main 2','Main 3','Main 4','Main 5','Post (adaptive)'},'FontSize',15)
ylabel('Accuracy (%)')
Name_legend = Names;
Name_legend{Nsub+1} = 'Mean';
legend(a,Name_legend,'Location','eastoutside')

%%
Nsub_sp = 8;
figure; hold on;
colors = gradient./255;
a = [];
for s= 1:Nsub_sp
a(s) = plot(Acc_online(s,:)','Color',colors(s,:),'marker','o','LineWidth',2,'MarkerFaceColor',colors(s,:));
end
a(Nsub_sp+1) = plot(mean(Acc_online,'omitnan')','Color','black','marker','o','LineWidth',4,'MarkerFaceColor','black');

plot([1.5 1.5],[0 105],'Color',[0.3 0.3 0.3],'LineStyle','--')
plot([6.5 6.5],[0 105],'Color',[0.3 0.3 0.3],'LineStyle','--')

xlim([0.5 7.5])
ylim([0 105])
set(gca,'XTick',1:7,'XTickLabel',{'Pre','Main 1','Main 2','Main 3','Main 4','Main 5','Post (adaptive)'},'FontSize',15)
ylabel('Accuracy (%)')
Name_legend = Names(1:Nsub_sp);
Name_legend{Nsub_sp+1} = 'Mean';
legend(a,Name_legend,'Location','eastoutside')
%---
Nsub_wat = 4;
figure; hold on;
colors = gradient./255;
a = [];
for s= 1:Nsub_wat
a(s) = plot(Acc_online(Nsub_sp+s,:)','Color',colors(Nsub_sp+s,:),'marker','o','LineWidth',2,'MarkerFaceColor',colors(Nsub_sp+s,:));
end
a(Nsub_wat+1) = plot(mean(Acc_online(Nsub_sp+1:end,:),'omitnan')','Color','black','marker','o','LineWidth',4,'MarkerFaceColor','black');

plot([1.5 1.5],[0 105],'Color',[0.3 0.3 0.3],'LineStyle','--')
plot([6.5 6.5],[0 105],'Color',[0.3 0.3 0.3],'LineStyle','--')

xlim([0.5 7.5])
ylim([0 105])
set(gca,'XTick',1:7,'XTickLabel',{'Pre','Main 1','Main 2','Main 3','Main 4','Main 5','Post (adaptive)'},'FontSize',15)
ylabel('Accuracy (%)')
Name_legend = Names(Nsub_sp+1:end);
Name_legend{Nsub_wat+1} = 'Mean';
legend(a,Name_legend,'Location','eastoutside')


%%
Acc_fixed_ = Acc_fixed;
Acc_fixed_(:,1) = Acc_online(:,1);
figure;
for s = 1:Nsub
    subplot(3,4,s);
    hold on;
    a1 = plot(Acc_online(s,:)','Color',colors(s,:),'marker','o','LineWidth',2,'MarkerFaceColor',colors(s,:));
    a2 = plot(Acc_fixed_(s,:)','Color',colors(s,:)*0.6,'marker','o','LineWidth',2,'MarkerFaceColor',colors(s,:),'LineStyle','--');
    plot([1.5 1.5],[0 105],'Color',[0.3 0.3 0.3],'LineStyle','--')
plot([6.5 6.5],[0 105],'Color',[0.3 0.3 0.3],'LineStyle','--')

xlim([0.5 7.5])
ylim([0 105])
set(gca,'XTick',1:7,'XTickLabel',{'Pre','Main 1','Main 2','Main 3','Main 4','Main 5','Post'},'FontSize',15)
ylabel('Accuracy (%)')
title(Names{s})
end
legend([a1 a2],{'Adaptive','Fixed'},'Location','eastoutside')

