close all; clear all

%%

LAMBDA1 = 0;
LAMBDA2 = 0;

%%

threshold1 = 0;
threshold2 = 0;
SS03_main2(LAMBDA1,LAMBDA2,threshold1,threshold2)


%%
u = load(['adaptiveSPHDCA\SS03_exp_adaptiveSPHDCA_LDA_PMean_NoThreNoLongVarNfVarUC_l1_',num2str(LAMBDA1),...
    '_l2_',num2str(LAMBDA2),'_thres1_',num2str(threshold1),'_thres2_',num2str(threshold2),'.mat']);

Acc_adMHDPA = squeeze(sum(u.outputs ==u.answers)./sum(~isnan(u.answers)))'*100;

%-- fixed 
[~,res2] = max(u.outputs2);
Acc2_uc = squeeze(sum(u.answers == squeeze(res2))./sum(~isnan(u.answers)))'*100;
%-- adaptive MHDPA1
[~,res1] = max(u.outputs1);
Acc1_uc = squeeze(sum(u.answers == squeeze(res1))./sum(~isnan(u.answers)))'*100;

%%
f = load(['adaptiveSPHDCA\SS03_exp_adaptiveSPHDCA_LDA_PMean_l1_0_l2_0_thres1_0.2_thres2_1.4.mat']);
[~,res2] = max(f.outputs2);
Acc2_f = squeeze(sum(f.answers == squeeze(res2))./sum(~isnan(f.answers)))'*100;

%%

% Define the start and end colors
start_color = [0, 0, 255]; % Blue
end_color = [255, 0, 0]; % Red

Names = {'Pilot11','Pilot12','Pilot13','Pilot14','Pilot16','Pilot17','Pilot18','Pilot19','Pilot20','Pilot21'};
Nsub = length(Names);
gradient = zeros(Nsub, 3); % Initialize the matrix for RGB values
for i = 1:3
    gradient(:, i) = linspace(start_color(i), end_color(i), Nsub);
end
colors = gradient./255;
Nsub = size(Acc1_uc,1);
figure;
for s = 1:Nsub
    subplot(3,4,s);
    hold on;
        a00 = plot(Acc2_f(s,:)','Color',colors(s,:)*0.5,'marker','^','LineWidth',2,'MarkerFaceColor','k','LineStyle','-.');    
    a0 = plot(Acc1_uc(s,:)','Color',colors(s,:)*0.5,'marker','o','LineWidth',2,'MarkerFaceColor',colors(s,:),'LineStyle',':');    
    a2 = plot(Acc2_uc(s,:)','Color',colors(s,:)*0.3,'marker','s','LineWidth',2,'MarkerFaceColor',colors(s,:));
    a3 = plot(Acc_adMHDPA(s,:)','Color',colors(s,:),'marker','diamond','LineWidth',2,'MarkerFaceColor',colors(s,:),'LineStyle','-');
    plot([1.5 1.5],[0 105],'Color',[0.3 0.3 0.3],'LineStyle','--')
plot([6.5 6.5],[0 105],'Color',[0.3 0.3 0.3],'LineStyle','--')

xlim([0.5 7.5])
ylim([0 105])
set(gca,'XTick',1:7,'XTickLabel',{'Pre','Main 1','Main 2','Main 3','Main 4','Main 5','Post'},'FontSize',15)
ylabel('Accuracy (%)')
title(Names{s})
end
legend([a00 a2 a0 a3],{'Fixed LW','Fixed MW','adaptive MHDPA1','adaptive MHDPA2'},'Location','eastoutside')
set(gcf,'Position',[ 329   283   1318   697])