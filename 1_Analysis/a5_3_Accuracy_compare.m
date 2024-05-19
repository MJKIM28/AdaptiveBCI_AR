%% Accuracy

SubNames = {'Subtest12','Subtest13','Subtest14','Subtest15','Subtest17',...
    'Subtest18','Subtest19','Subtest20','Subtest21','Subtest22'};
Names = {'Pilot11','Pilot12','Pilot13','Pilot14','Pilot16','Pilot17','Pilot18','Pilot19','Pilot20','Pilot21'};


%% Fixed 
Acc_fix = [0.8	0.8	0.6667	0.8	0.8667	0.8	0.8667;
1	0.6667	0.7333	0.7333	0.7333	0.8	0.8667;
0.8	0.1429	0.0667	0.2	0.4	0.6	0.6667;
0.6	0.6667	0.5333	0.8667	0.6	0.8	0.8;
 1	NaN	1	1	0.9286	0.9286	1;
0.9286	0.8571	0.9286	0.5714	0.9286	0.8571	1;
0.7143	0.5714	0.8571	0.9286	0.6429	0.8571 NaN;	
0.9286	0.6429	0.6429	0.5714	0.7143	0.6429 NaN;
0.8571	0.3571	0.5714	0.3571	0.5714	0.2857	0.8571;
1	0.7857	0.9286	0.6429	0.7143	0.7857	1]*100;

%% Adaptive: select trials
Acc_ad_seltr = [0.8	0.7333	0.5333	0.7333	0.6667	0.6667		0.8667;
1	0.6	0.7333	0.6667	0.8	0.6667		0.8667;
0.8	0.1429	0.1333	0.2667	0.4	0.5333		0.5333;
0.6	0.6667	0.4667	0.7333	0.7333	0.7333		0.7333;
1 NaN	1	1	0.9286	1	1;
0.9286	0.8571	0.9286	0.6429	0.7857	0.7857	0.9286;
0.7143	0.5	1	0.8571	0.7143	0.7143	NaN;
0.9286	0.7143	0.6429	0.5	0.6429	0.6429 NaN	;
0.8571	0.2857	0.5	0.2143	0.4286	0.2857	0.5714;
1	0.8571	0.7857	0.7857	0.7857	0.8571	0.9286]*100;

%% Adaptive: select blocks (remove  oldest nonSV)
Acc_ad_selbl = [0.8	0.8	0.7333	0.8	0.7333	0.7333	0.8667;
1	0.6667	0.8	0.8	0.6	0.6667	0.6667;
0.8	0.1429	0.2	0.2667	0.3333	0.3333	0.6;
0.6	0.6667	0.4667	0.9333	0.6667	0.6	0.9333;
1	NaN	1	1	1	1		1;
0.9286	0.7857	0.8571	0.5714	0.7143	0.6429		1;
0.7143	0.5714	0.9286	0.8571	0.7857	0.7143	NaN;	
0.9286	0.7143	0.6429	0.6429	0.7143	0.6429	NaN;	
0.8571	0.4286	0.3571	0.4286	0.4286	0.3571		0.7143;
1	0.8571	0.8571	0.6429	0.6429	0.6429		0.7857]*100;

%% Adaptive: select blocks (remove oldest)
Acc_ad_selbl_remold = [0.8	0.9333    0.9333    0.7333    0.7333    0.6667    0.8667;
1	0.8000    0.6000    0.9333    0.6667    0.5333    0.7333;
0.8	 0.1429    0.3333    0.2000    0.6000    0.3333    0.6667;
0.6	  0.8000    0.5333    0.8000    0.7333    0.7333    0.8000;
1	NaN	 1.0000    0.9286    1.0000    1.0000    1.0000;
0.9286	0.7857    0.8571    0.6429    0.8571    0.7857    1.0000;
0.7143	0.6429    1.0000    0.9286    0.8571    0.8571	NaN;	
0.9286	0.7857    0.7143    0.7143    0.7857    0.8571	NaN;	
0.8571	0.3571    0.3571    0.4286    0.5000    0.3571    0.9286;
1	0.8571    0.9286    0.7857    0.5000    0.6429    1.0000]*100;

save('Accuracy_BCI_exp','Acc_fix','Acc_ad_seltr','Acc_ad_selbl');
%% adpative MHDPA
lambda1 = 0.03; lamba2 = 0.04;
thres1 = 0.3; thres2 = 0.3;
load(['adaptiveSPHDCA\SS03_exp_adaptiveSPHDCA_LDA_l1_0.03_l2_0.03_thres1_0.3_thres2_0.3.mat'])
Acc_adMHDPA = squeeze(sum(outputs ==answers)./sum(~isnan(answers)))'*100;


%% adaptation methods
Nsub = length(Names);

% Define the start and end colors
start_color = [0, 0, 255]; % Blue
end_color = [255, 0, 0]; % Red

gradient = zeros(Nsub, 3); % Initialize the matrix for RGB values
for i = 1:3
    gradient(:, i) = linspace(start_color(i), end_color(i), Nsub);
end
colors = gradient./255;

figure;
for s = 1:Nsub
    subplot(3,4,s);
    hold on;
    a0 = plot(Acc_ad_selbl(s,:)','Color',colors(s,:),'marker','o','LineWidth',2,'MarkerFaceColor',colors(s,:));    
    a1 = plot(Acc_ad_seltr(s,:)','Color',colors(s,:),'marker','^','LineWidth',2,'MarkerFaceColor',colors(s,:),'LineStyle',':');
    a2 = plot(Acc_fix(s,:)','Color',colors(s,:)*0.3,'marker','s','LineWidth',2,'MarkerFaceColor',colors(s,:));
     a3 = plot(Acc_adMHDPA(s,:)','Color',colors(s,:)*0.5,'marker','diamond','LineWidth',2,'MarkerFaceColor',colors(s,:),'LineStyle','--');
    plot([1.5 1.5],[0 105],'Color',[0.3 0.3 0.3],'LineStyle','--')
plot([6.5 6.5],[0 105],'Color',[0.3 0.3 0.3],'LineStyle','--')

xlim([0.5 7.5])
ylim([0 105])
set(gca,'XTick',1:7,'XTickLabel',{'Pre','Main 1','Main 2','Main 3','Main 4','Main 5','Post'},'FontSize',15)
ylabel('Accuracy (%)')
title(Names{s})
end
legend([a0 a1 a3 a2],{'Ad(block)','Ad(trial)','AdMHDPA','Fixed'},'Location','eastoutside')
set(gcf,'Position',[ 329   283   1318   697])
%% MHDPA
%-- fixed LDA
[~,res2] = max(outputs2);
Acc2 = squeeze(sum(answers == squeeze(res2))./sum(~isnan(answers)))'*100;
%-- adaptive MHDPA1
[~,res1] = max(outputs1);
Acc1 = squeeze(sum(answers == squeeze(res1))./sum(~isnan(answers)))'*100;


figure;
for s = 1:Nsub
    subplot(3,4,s);
    hold on;
    a0 = plot(Acc1(s,:)','Color',colors(s,:)*0.5,'marker','o','LineWidth',2,'MarkerFaceColor',colors(s,:),'LineStyle',':');    
    a2 = plot(Acc2(s,:)','Color',colors(s,:)*0.3,'marker','s','LineWidth',2,'MarkerFaceColor',colors(s,:));
    a3 = plot(Acc_adMHDPA(s,:)','Color',colors(s,:),'marker','diamond','LineWidth',2,'MarkerFaceColor',colors(s,:),'LineStyle','-');
    plot([1.5 1.5],[0 105],'Color',[0.3 0.3 0.3],'LineStyle','--')
plot([6.5 6.5],[0 105],'Color',[0.3 0.3 0.3],'LineStyle','--')

xlim([0.5 7.5])
ylim([0 105])
set(gca,'XTick',1:7,'XTickLabel',{'Pre','Main 1','Main 2','Main 3','Main 4','Main 5','Post'},'FontSize',15)
ylabel('Accuracy (%)')
title(Names{s})
end
legend([a2 a0 a3],{'Fixed DSP-LDA','adaptive MHDPA1','adaptive MHDPA2'},'Location','eastoutside')
set(gcf,'Position',[ 329   283   1318   697])

%% Distraction effect
AccDiff1 = Acc_ad_seltr(1:4,2:end-1) - Acc_ad_seltr(1:4,1)

AccDiff2 = Acc_ad_selbl(5:end,2:end-1) - Acc_ad_selbl(5:end,1)

Xs = [mean(AccDiff1,2,'omitnan'); mean(AccDiff2,2,'omitnan')]

[p,h,stat] = signrank(Xs)


figure; 
violinplot(Xs)
xlim([ 0 2])
set(gca,'xticklabel','','fontsize',15)
ylabel('\Delta Accuracy')