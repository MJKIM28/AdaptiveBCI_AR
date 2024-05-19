close all; clear all

%%

LAMBDA1 = 0;
LAMBDA2 = 0;

%%

threshold1 = 0.3;
threshold2 = 0.3;
SS03_main(LAMBDA1,LAMBDA2,threshold1,threshold2)

%%
uc_gc = load(['adaptiveSPHDCA\SS03_exp_adaptiveSPHDCA_LDA_PMean_GC_l1_',num2str(LAMBDA1),...
    '_l2_',num2str(LAMBDA2),'_thres1_',num2str(threshold1),'_thres2_',num2str(threshold2),'.mat']);

Acc_adMHDPA_ucgc = squeeze(sum(uc_gc.outputs ==uc_gc.answers)./sum(~isnan(uc_gc.answers)))'*100;

%%
uc_ft = load(['adaptiveSPHDCA\SS03_exp_adaptiveSPHDCA_LDA_PMean_l1_',num2str(LAMBDA1),...
    '_l2_',num2str(LAMBDA2),'_thres1_',num2str(threshold1),'_thres2_',num2str(threshold2),'.mat']);

Acc_adMHDPA_uc = squeeze(sum(uc_ft.outputs ==uc_ft.answers)./sum(~isnan(uc_ft.answers)))'*100;
%-- fixed LDA
[~,res2] = max(uc_ft.outputs2);
Acc2_uc = squeeze(sum(uc_ft.answers == squeeze(res2))./sum(~isnan(uc_ft.answers)))'*100;
%-- adaptive MHDPA1
[~,res1] = max(uc_ft.outputs1);
Acc1_uc = squeeze(sum(uc_ft.answers == squeeze(res1))./sum(~isnan(uc_ft.answers)))'*100;

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
    a0 = plot(Acc1_uc(s,:)','Color',colors(s,:)*0.5,'marker','o','LineWidth',2,'MarkerFaceColor',colors(s,:),'LineStyle',':');    
    a2 = plot(Acc2_uc(s,:)','Color',colors(s,:)*0.3,'marker','s','LineWidth',2,'MarkerFaceColor',colors(s,:));
    a3 = plot(Acc_adMHDPA_uc(s,:)','Color',colors(s,:),'marker','diamond','LineWidth',2,'MarkerFaceColor',colors(s,:),'LineStyle','-');
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
%% update coefficient function
figure;
sigma = 0.2;
A = 0.1;
mu = 0;
x = -1:0.01:1;
f = A * exp(-(x - mu).^2 / (2 * sigma^2));
plot(x, f, '-r', 'LineWidth', 2);
xlabel('Cosine similarity')
ylabel('Update coefficient')
set(gca,'fontsize',20)
title('DSP update coefficient rule')

figure;
A = 0.06;
f = A * exp(-(x - mu).^2 / (2 * sigma^2));
plot(x, f, '-r', 'LineWidth', 2);
xlabel('Cosine similarity')
ylabel('Update coefficient')
set(gca,'fontsize',20)
title('LDA update coefficient rule')
%%
fc = load(['adaptiveSPHDCA\SS03_exp_adaptiveSPHDCA_LDA_l1_0.03_l2_0.03_thres1_0.3_thres2_0.3.mat']);
%%
Acc_adMHDPA = squeeze(sum(fc.outputs ==fc.answers)./sum(~isnan(fc.answers)))'*100;
%-- fixed LDA
[~,res2] = max(fc.outputs2);
Acc2 = squeeze(sum(fc.answers == squeeze(res2))./sum(~isnan(fc.answers)))'*100;
%-- adaptive MHDPA1
[~,res1] = max(fc.outputs1);
Acc1 = squeeze(sum(fc.answers == squeeze(res1))./sum(~isnan(fc.answers)))'*100;



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

%%
colorcode = [0 0 0;0, 0, 128;107, 142, 35; 255, 20, 147]./255;

figure;
for s = 1:Nsub
    subplot(3,4,s);
    hold on;
    a1 = plot(Acc2(s,:)','Color','k','marker','s','LineWidth',2,'MarkerFaceColor','k','LineStyle',':');   
    a2 = plot(Acc_adMHDPA(s,:)','Color',colorcode(2,:),'marker','o','LineWidth',2,'MarkerFaceColor',colorcode(2,:),'LineStyle','-');
     a3 = plot(Acc_adMHDPA_ucgc(s,:)','Color',colorcode(3,:),'marker','^','LineWidth',2,'MarkerFaceColor',colorcode(3,:),'LineStyle','-');  
    a4 = plot(Acc_adMHDPA_uc(s,:)','Color',colorcode(4,:),'marker','diamond','LineWidth',2,'MarkerFaceColor',colorcode(4,:),'LineStyle','-');
    
    plot([1.5 1.5],[0 105],'Color',[0.3 0.3 0.3],'LineStyle','--')
plot([6.5 6.5],[0 105],'Color',[0.3 0.3 0.3],'LineStyle','--')

xlim([0.5 7.5])
ylim([0 105])
set(gca,'XTick',1:7,'XTickLabel',{'Pre','Main 1','Main 2','Main 3','Main 4','Main 5','Post'},'FontSize',15)
ylabel('Accuracy (%)')
title(Names{s})
end
legend([a1 a2 a3 a4],{'Fixed DSP-LDA','adaptive MHDPA (PM,GC)','adaptive MHDPA (PM,GC, w/ aUC)','adaptive MHDPA (PM, w/ aUC)'},'Location','eastoutside')
set(gcf,'Position',[ 329   283   1318   697])
