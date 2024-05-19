
LAMBDA1 = 0;
LAMBDA2 = 0;

threshold1 = 0;
threshold2 = 0;

SS03_main(LAMBDA1,LAMBDA2,threshold1,threshold2)

%%
uc = load(['adaptiveSPHDCA\SS03_exp_adaptiveSPHDCA_LDA_PMean_BaseThres_l1_',num2str(LAMBDA1),...
    '_l2_',num2str(LAMBDA2),'_thres1_',num2str(threshold1),'_thres2_',num2str(threshold2),'.mat']);
%