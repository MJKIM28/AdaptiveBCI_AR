
%% Accuracy - train again
%%
for s = 1:Nsub
    SubName = SubNameList{s};

    load([fpath,'\',SubName,'.mat'])
    targetlist = vars.target;
    %% load accuracy
    p = load([fpath,'\Dat_',SubName,'\param.mat']);

    pred_pre =  p.param.prediction(end-Ntr_te+1:end-(Nsess-1)*Ntr_con);
    pred_main = p.param.prediction(end-(Nsess-1)*Ntr_con+1:end-(Nsess-5)*Ntr_con);
    pred_post = p.param.prediction(end-(Nsess-5)*Ntr_con+1:end);

    vars.Acc_pre = mean(pred_pre' == targetlist(:,1));

    target_main =  targetlist(:,2:5);
    vars.Acc_mainAll = mean(pred_main'==target_main(:));

    target_post = targetlist(:,6:end);
    vars.Acc_post = mean(pred_post'==target_post(:));

    for s= 1:Nsess-1-1
        vars.Acc_main_sess(s) = mean(pred_main(Ntr_con*(s-1)+1:s*Ntr_con)' == targetlist(:,1+s));
    end
    fprintf('==============\nAcc_pre  %.2f \n \nAcc_main1  %.2f \nAcc_main2  %.2f\nAcc_main3  %.2f\nAcc_main4  %.2f\nAcc_main(all) %.2f\n\nAcc_post %.2f\n',...
        vars.Acc_pre,vars.Acc_main_sess(1),vars.Acc_main_sess(2),vars.Acc_main_sess(3),vars.Acc_main_sess(4),vars.Acc_mainAll,vars.Acc_post);

    Acc_online = [vars.Acc_pre vars.Acc_main_sess vars.Acc_post];

    %% re-estimate (fixed)
    %-- train
    cd(codepath)
    targets = vars.target(:);
%% Performance re-generation
% 1. train again
% 2. use classifier from exp

clear all; close all
%% parameters
fpath = 'E:\1_연구\분석코드\2024_adaptive_BCI\Data\Pilot';
codepath = 'E:\1_연구\분석코드\2024_adaptive_BCI\AdaptiveBCI_AR\0_online_code';

Ntr_tr = 8;
Ntr_con = 15;
Nsess = 6;
Ntr_te = Ntr_con*Nsess;

B_fir_lf_2 = getfirfiltcoeff(15,'low',500,0);
LPF15 = 1;

% SubNameList = {'Subtest02','Subtest03','Subtest04','Subtest05'};
SubNameList = {'Subtest06','Subtest07'};

Nsub = length(SubNameList);
    e = load([fpath,'\Train\Epoch\',SubName]);

    param = p.param;
    param.trD.mode = 'training';
    [Feature, label, param] = FeatureExtraction(e.Epoch, param);
    [C,param_in]               = Classification(Feature,label, param);


    %-- test
    et = load([fpath,'\Test\Epoch\',SubName]);
    pred = [];
    for tr = 1:Ntr_te
        [Feature, label, param_in] = FeatureExtraction(et.Epoch{tr}, param_in);
        [C,param_in]               = Classification(Feature,label, param_in);
        pred(tr) = C;
    end

    targets_con = reshape(targets,[Ntr_con,Nsess]);
    pred_con = reshape(pred,[Ntr_con,Nsess]);

    Acc_fixed = [];
    for con = 1:Nsess
        Acc                                                              _fixed(con) = mean(targets_con(:,con) == pred_con(:,con));
    end

    %% re-estimate (adpative)
    % pred_ad = [];
    % for tr = 1:Ntr_te
    %
    %
    %     paramNew = p.param;
    %     if tr > 1
    %         if isfield(p.param.update{tr-1},'DSP')
    %             paramNew.DSP = p.param.update{tr-1}.DSP;
    %             paramNew.trD.mdl = p.param.update{tr-1}.mdl.mdl;
    %         else
    %             paramNew.DSP = p.param.update{tr-2}.DSP;
    %             paramNew.trD.mdl = p.param.update{tr-2}.mdl.mdl;
    %         end
    %     end
    %     [Feature, label, paramNew] = FeatureExtraction(et.Epoch{tr}, paramNew);
    %     C               = Classification(Feature,label, paramNew);
    %     pred_ad(tr) = C;
    %
    % end
    %
    % pred_ad_con = reshape(pred_ad,[Ntr_con,Nsess]);
    %
    % Acc_ad = [];
    % for con = 1:Nsess
    %     Acc_ad(con) = mean(targets_con(:,con) == pred_ad_con(:,con));
    % end

    mkdir([fpath,'\Accuracy'])
    save([fpath,'\Accuracy\',SubName],'Acc_fixed','Acc_online');

end
%% Accuracy - use classifier from exp

