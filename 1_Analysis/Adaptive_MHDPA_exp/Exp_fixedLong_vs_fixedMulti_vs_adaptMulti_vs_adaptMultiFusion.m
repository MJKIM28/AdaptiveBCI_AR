%%

% 1. fixed long win
LDA_LongWin_main('fixed','NotFusion')
% 2. fixed multi
% LDA_multiOnly_main('fixed','NotFusion')

% 3. adaptive multi
LDA_multiOnly_main('adaptive','NotFusion')

% 4. adaptive-fixed fusion multi
LDA_multiOnly_main('adaptive','fusion')

%%

f_l = load('adaptiveMWHDPA_LDA_LongWin_NoThre_VarNf_VarUC_original\fixed_NotFusion.mat');
[~,res2] = max(f_l.outputs2);
Acc_f_l = squeeze(sum(squeeze(res2) ==f_l.answers)./sum(~isnan(f_l.answers)))'*100; % fixed long
[~,res1] = max(f_l.outputs1);
Acc_f_m = squeeze(sum(squeeze(res1) ==f_l.answers)./sum(~isnan(f_l.answers)))'*100; % fixed multi


a_nf = load('adaptiveMWHDPA_LDA_MultiONly_NoThre_VarNf_VarUC_original\adaptive_NotFusion.mat');
[~,res1] = max(a_nf.outputs1);
Acc_a_nf = squeeze(sum(squeeze(res1) ==a_nf.answers)./sum(~isnan(a_nf.answers)))'*100; % adapt not fusion

a_f = load('adaptiveMWHDPA_LDA_MultiONly_NoThre_VarNf_VarUC_original\adaptive_fusion.mat');
Acc_a_f = squeeze(sum(a_f.outputs ==a_f.answers)./sum(~isnan(a_f.answers)))'*100; % adapt fusion



%%


Names = {'Pilot11','Pilot12','Pilot13','Pilot14','Pilot16','Pilot17','Pilot18','Pilot19','Pilot20','Pilot21'};
Nsub = length(Names);
colors = [125 125 125;0, 0, 128;107, 142, 35; 255, 20, 147]./255;

figure;
for s = 1:Nsub
    subplot(3,4,s);
    hold on;
        a00 = plot(Acc_f_l(s,:)','Color',colors(1,:),'marker','^','LineWidth',2,'MarkerFaceColor',colors(1,:),'LineStyle','-.');    
    a0 = plot(Acc_f_m(s,:)','Color',colors(2,:),'marker','o','LineWidth',2,'MarkerFaceColor',colors(2,:),'LineStyle',':');    
    a2 = plot(Acc_a_nf(s,:)','Color',colors(3,:),'marker','s','LineWidth',2,'MarkerFaceColor',colors(3,:));
    a3 = plot(Acc_a_f(s,:)','Color',colors(4,:),'marker','diamond','LineWidth',2,'MarkerFaceColor',colors(4,:),'LineStyle','-');
    plot([1.5 1.5],[0 105],'Color',[0.3 0.3 0.3],'LineStyle','--')
plot([6.5 6.5],[0 105],'Color',[0.3 0.3 0.3],'LineStyle','--')

xlim([0.5 7.5])
ylim([0 105])
set(gca,'XTick',1:7,'XTickLabel',{'Pre','Main 1','Main 2','Main 3','Main 4','Main 5','Post'},'FontSize',15)
ylabel('Accuracy (%)')
title(Names{s})
end
legend([a00 a0 a2 a3],{'Fixed LW','Fixed MW','adaptive MW not fusion','adaptive MW fusion'},'Location','eastoutside')
set(gcf,'Position',[ 329   283   1318   697])


%% update LR until 400sample (10 block)

% 3. adaptive multi
LDA_multiOnly_main('adaptive','NotFusion')

% 4. adaptive-fixed fusion multi
LDA_multiOnly_main('adaptive','fusion')

%%
a_nf_ulr = load('adaptiveMWHDPA_LDA_MultiONly_NoThre_VarNf_VarUC_updateLRmax50\adaptive_NotFusion.mat');
[~,res1] = max(a_nf_ulr.outputs1);
Acc_a_nf_ulr = squeeze(sum(squeeze(res1) ==a_nf_ulr.answers)./sum(~isnan(a_nf_ulr.answers)))'*100; % adapt not fusion

a_f_ulr = load('adaptiveMWHDPA_LDA_MultiONly_NoThre_VarNf_VarUC_updateLRmax50\adaptive_fusion.mat');
Acc_a_f_ulr = squeeze(sum(a_f_ulr.outputs ==a_f_ulr.answers)./sum(~isnan(a_f_ulr.answers)))'*100; % adapt fusion


Names = {'Pilot11','Pilot12','Pilot13','Pilot14','Pilot16','Pilot17','Pilot18','Pilot19','Pilot20','Pilot21'};
Nsub = length(Names);
colors = [125 125 125;0, 0, 128;107, 142, 35; 255, 20, 147]./255;

figure;
for s = 1:Nsub
    subplot(3,4,s);
    hold on;
        a00 = plot(Acc_f_l(s,:)','Color',colors(1,:),'marker','^','LineWidth',2,'MarkerFaceColor',colors(1,:),'LineStyle','-.');    
    a0 = plot(Acc_f_m(s,:)','Color',colors(2,:),'marker','o','LineWidth',2,'MarkerFaceColor',colors(2,:),'LineStyle',':');    
    a2 = plot(Acc_a_nf_ulr(s,:)','Color',colors(3,:),'marker','s','LineWidth',2,'MarkerFaceColor',colors(3,:));
    a3 = plot(Acc_a_f(s,:)','Color',colors(4,:),'marker','diamond','LineWidth',2,'MarkerFaceColor',colors(4,:),'LineStyle','-');
    plot([1.5 1.5],[0 105],'Color',[0.3 0.3 0.3],'LineStyle','--')
plot([6.5 6.5],[0 105],'Color',[0.3 0.3 0.3],'LineStyle','--')

xlim([0.5 7.5])
ylim([0 105])
set(gca,'XTick',1:7,'XTickLabel',{'Pre','Main 1','Main 2','Main 3','Main 4','Main 5','Post'},'FontSize',15)
ylabel('Accuracy (%)')
title(Names{s})
end
legend([a00 a0 a2 a3],{'Fixed LW','Fixed MW','adaptive MW not fusion','adaptive MW fusion'},'Location','eastoutside')
set(gcf,'Position',[ 329   283   1318   697])
