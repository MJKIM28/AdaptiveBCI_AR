clear all; close all;
%% Accuracy - train again

fpath = 'E:\1_연구\분석코드\2024_adaptive_BCI\Data\Pilot';
codepath = 'E:\1_연구\분석코드\2024_adaptive_BCI\AdaptiveBCI_AR\0_online_code';


SubNameList = {'Subtest02','Subtest03','Subtest04','Subtest05','Subtest06','Subtest07'};

Nsub = length(SubNameList);

Ntr_tr = 8;
Ntr_con = 15;
Nsess = 6;
Ntr_te = Ntr_con*Nsess;


%%
for s = setdiff(1:Nsub,5)
    SubName = SubNameList{s};

    load([fpath,'\',SubName,'.mat'])
    targetlist = vars.target;
    targetlist(targetlist == 0) = NaN;

    trials = find(~isnan(targetlist(:)));
    targets = targetlist;

    %% load accuracy
%     p = load([fpath,'\Dat_',SubName,'\param.mat']);
% 
%     pred_pre =  p.param.prediction(end-Ntr_te+1:end-(Nsess-1)*Ntr_con);
%     pred_main = p.param.prediction(end-(Nsess-1)*Ntr_con+1:end-(Nsess-5)*Ntr_con);
%     pred_post = p.param.prediction(end-(Nsess-5)*Ntr_con+1:end);
% 
%     vars.Acc_pre = mean(pred_pre' == targetlist(:,1));
% 
%     target_main =  targetlist(:,2:5);
%     vars.Acc_mainAll = mean(pred_main'==target_main(:));
% 
%     target_post = targetlist(:,6:end);
%     vars.Acc_post = mean(pred_post'==target_post(:));
% 
%     for s= 1:Nsess-1-1
%         vars.Acc_main_sess(s) = mean(pred_main(Ntr_con*(s-1)+1:s*Ntr_con)' == targetlist(:,1+s));
%     end
%     fprintf('==============\nAcc_pre  %.2f \n \nAcc_main1  %.2f \nAcc_main2  %.2f\nAcc_main3  %.2f\nAcc_main4  %.2f\nAcc_main(all) %.2f\n\nAcc_post %.2f\n',...
%         vars.Acc_pre,vars.Acc_main_sess(1),vars.Acc_main_sess(2),vars.Acc_main_sess(3),vars.Acc_main_sess(4),vars.Acc_mainAll,vars.Acc_post);
% 
%     Acc_online = [vars.Acc_pre vars.Acc_main_sess vars.Acc_post];

    %% re-estimate (fixed)
    %-- train
    cd(codepath)
%% Performance re-generation
% 1. train again
% 2. use classifier from exp

 close all
%% parameters


    e = load([fpath,'\Train\Epoch\',SubName]);
    p = load([fpath,'\Dat_',SubName,'\param.mat']);
    param = p.param;

    if isfield(param.trD,'mdl_init') && isfield(param.DSP,'init')
        fprintf(['%s : use original model\n'],SubName);
        param.trD.mdl = param.trD.mdl_init;
        param.trD.mdl_adapt = param.trD.mdl; % initialization

        lambda = param.DSP.lambda;
        param.DSP = param.DSP.init;
        param.DSP.init = param.DSP;
        param.DSP.lambda = lambda;

        param.trD.feature = param.trD.feature(1:320,:);
        param.trD.label = param.trD.label(1:320,:);
        
        param_in = param;

    else
        fprintf(['%s : train from scratch\n'],SubName);

        param.trD.mode = 'training';
    [Feature, label, param] = FeatureExtraction(e.Epoch, param);
    [C,param_in]               = Classification(Feature,label, param);
    param_in.trD.mdl_init = param_in.trD.mdl;
    param_in.trD.mdl_adapt = param_in.trD.mdl;
    param_in.trD.feature = Feature;
    param_in.trD.label = label;

    end

    X_init{1} = e.Epoch.tar;
    X_init{2} = e.Epoch.nar;
        
    


    %-- pre
    et = load([fpath,'\Test\Epoch\',SubName]);
    pred = [];
    N_te = length(et.Epoch);
    for tr = 1:N_te 
        [Feature, label, param_in] = FeatureExtraction(et.Epoch{tr}, param_in);
        [C,param_in]               = Classification(Feature,label, param_in);
        pred(tr) = C;
    end

    targets_con = reshape(targets,Ntr_con,[]);
    pred_con = reshape(pred,Ntr_con,[]);

    
    Acc_fixed= mean(targets_con== pred_con);
    

    %% re-estimate (adpative)
    tt = 1;

    UPmodel = [];
    pred_ad = [];
    param_ad = param_in;
    X = X_init;
    for tr = Ntr_con+1:N_te
    
    %-- use saved model
%         paramNew = p.param;
%         if tr > 1
%             if isfield(p.param.update{tr-1},'DSP')
%                 paramNew.DSP = p.param.update{tr-1}.DSP;
%                 paramNew.trD.mdl = p.param.update{tr-1}.mdl.mdl;
%             else
%                 paramNew.DSP = p.param.update{tr-2}.DSP;
%                 paramNew.trD.mdl = p.param.update{tr-2}.mdl.mdl;
%             end
%         end

%-- re-train

        param_ad.trD.mode = 'testing';
        param_ad.trD.mdl = param_ad.trD.mdl_adapt;
        [Feature, label, param_ad] = FeatureExtraction(et.Epoch{tr}, param_ad);
        [C,param_ad] = Classification(Feature,label, param_ad);

      [upcheck,up,Posterior,trs,D] = checkupdate(param_ad.trD.score, [] ,param_ad.trD.threshold,C,[]);


      UPmodel(tt).upcheck = upcheck;
      UPmodel(tt).upids = up;
      UPmodel(tt).Posterior = Posterior;
      UPmodel(tt).trials_candidate = trs;
      UPmodel(tt).Dist = D;


      if upcheck && ~isempty(up)
          EP_1block = Epoch_condition(et.Epoch{tr},param_ad);

          %-- update DSP
          EP_update = EP_1block;
          EP_update.nar = EP_update.nar(:,:,up,:);

          %-- update DSP
%           param_ad  = updateDSP(EP_update,C,param_ad);

          %-- apply updated DSP to feat

          Repeat = param_ad.repeat;
          param_ad.repeat = size(EP_update.nar,3);
          Feat = FeatureExt_DSP(EP_update,param_ad);
          [C_new,param_ad] = Classification(Feat,[],param_ad);
% 
%           label_temp = -ones(size(Feat,1),1);
%           label_temp = reshape(label_temp,param_ad.repeat,param_ad.NumStims);
%           label_temp(:,C_new) = 1;
%           label = label_temp(:);
% 
%           %-- re-calibration
%           param_ad.trD.feature = [param_ad.trD.feature; Feat];
%           param_ad.trD.label = [param_ad.trD.label; label];
%           [~,param_ad] = Classification(param_ad.trD.feature,param_ad.trD.label,param_ad);

%           param_ad.repeat = Repeat;
          %--- 
          %-- re-estimate DSP & SVM
          
        erptar = EP_update.nar(:,:,:,C);
        erpnar = EP_update.nar(:,:,:,setdiff(1:param.NumStims,C));
        X{1} = cat(3,X{1},erptar);
        X{2} = cat(3,X{2},erpnar(:,:,:));
        EP_new.tar = X{1};
        EP_new.nar = X{2};
        param_ad.trD.mode = 'training';
        param_ad.DSP = rmfield(param_ad.DSP,'W');
        [feat_new,label_new, param_ad] = FeatureExt_DSP(EP_new, param_ad);
        [~,param_ad] = Classification(feat_new,label_new,param_ad);

        %-- 


          fprintf('>> Updated\n')


          UPmodel(tt).DSP = param_ad.DSP;
          UPmodel(tt).mdl = param_ad.trD;
          param_ad.trD.mdl_adapt = param_ad.trD.mdl;
      else
          fprintf('>> Not updated\n')
      end


        pred_ad(tt) = C;
    
        tt = tt + 1;
    end
    
    pred_ad_con = reshape(pred_ad,Ntr_con,[]);

    Acc_ad = mean(targets_con(:,2:end) == pred_ad_con);
    

    mkdir([fpath,'\Simulate_reestDSP_Acc'])
    save([fpath,'\Simulate_reestDSP_Acc\',SubName],'Acc_fixed','Acc_ad');

    mkdir([fpath,'\Simulate_reestDSP_param'])
    save([fpath,'\Simulate_reestDSP_param\',SubName],'UPmodel');

end
%%
Acc_fix_all = []; Acc_ad_all = [];
for s = setdiff(1:Nsub,5)
        SubName = SubNameList{s};
    load([fpath,'\Simulate_reestDSP_Acc\',SubName])
    nsess = size(Acc_fixed,2);
 Acc_fix_all(s,1:nsess) = Acc_fixed;
 Acc_ad_all(s,1:nsess-1) = Acc_ad;
end