clear all; close all;
%% Accuracy - train again

fpath = 'E:\1_연구\분석코드\2024_adaptive_BCI\Data\Pilot';
codepath = 'E:\1_연구\분석코드\2024_adaptive_BCI\AdaptiveBCI_AR\0_online_code';


SubNameList = {'Subtest02','Subtest03','Subtest04','Subtest05','Subtest07'};

Nsub = length(SubNameList);

Ntr_tr = 8;
Ntr_con = 15;
Nsess = 6;
Ntr_te = Ntr_con*Nsess;

LAMBDA1 = 0.01;
LAMBDA2 = 0.04;
REMOVE = 1;
%%
for s = 1:Nsub
    SubName = SubNameList{s};

    load([fpath,'\',SubName,'.mat'])
    targetlist = vars.target;
    targetlist(targetlist == 0) = NaN;

    trials = find(~isnan(targetlist(:)));
    targets = targetlist(:);


    %% re-estimate (fixed)
    %-- train
    cd(codepath)
%% Performance re-generation
% 1. train again
% 2. use classifier from exp

 close all
%% parameters
rng(1)

    e = load([fpath,'\Train\Epoch\',SubName]);
    p = load([fpath,'\Dat_',SubName,'\param.mat']);
    param = p.param;

    
        fprintf(['%s : train from scratch\n'],SubName);

        param.trD.mode = 'training';
         EP = Epoch_condition(e.Epoch,param);
         param1 = param;
         param2 = param;
         param1.DSP.window = param.Baseline + 1:5:param.Baseline+0.3*param.Fs;
         param2.DSP.window = param.Baseline + 0.3*param.Fs+1:5:param.Baseline + param.Epocline;
         param1.DSP.nf = 3; param2.DSP.nf = 3;
         param1.DSP = rmfield(param1.DSP,'W');
         param2.DSP = rmfield(param2.DSP,'W');         
         [feat1,label,param1] = FeatureExt_DSP(EP,param1);  
         [feat2,~,param2] = FeatureExt_DSP(EP,param2);    
         Feature = [feat1 feat2];

         [C,param_in]               = Classification(Feature,label, param);
    
         param_in.trD.mdl_init = param_in.trD.mdl;
    param_in.trD.mdl_adapt = param_in.trD.mdl;
    param_in.trD.feature = Feature;
    param_in.trD.label = label;

    

    X_init{1} = e.Epoch.tar;
    X_init{2} = e.Epoch.nar;
        
    

    param1.trD.mode = 'testing';
    param2.trD.mode = 'testing';
    %-- pre
    et = load([fpath,'\Test\Epoch\',SubName]);
    pred = [];
    N_te = length(et.Epoch);
    for tr = 1:N_te 
                 EP = Epoch_condition(et.Epoch{tr},param_in);
                 [feat1,label,param1] = FeatureExt_DSP(EP,param1);  
         [feat2,~,param2] = FeatureExt_DSP(EP,param2);    
         Feature = [feat1 feat2];

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
    param1_ad = param1;
    param2_ad = param2;
    param_ad.trD.threshold = 0.4;
    X = X_init;
%     UPmodel(1).DSP_init = param_in.DSP;
%     UPmodel(1).mdl_init = param_in.trD.mdl;
    for tr = Ntr_con+1:N_te
    rng(1)


%-- re-train
param1_ad.DSP.lambda = LAMBDA1;
param2_ad.DSP.lambda = LAMBDA2;

        param_ad.trD.mode = 'testing';
        param_ad.trD.mdl = param_ad.trD.mdl_adapt;


        EP = Epoch_condition(et.Epoch{tr},param_ad);
        [feat1,label,param1] = FeatureExt_DSP(EP,param1_ad);
        [feat2,~,param2] = FeatureExt_DSP(EP,param2_ad);
        Feature = [feat1 feat2];

        [C,param_ad]               = Classification(Feature,label, param_ad);



      [upcheck,up,Posterior,trs,D] = checkupdate(param_ad.trD.score, [] ,param_ad.trD.threshold,C,[]);


%       UPmodel(tt).classprob = param_ad.trD.score;
%       UPmodel(tt).upcheck = upcheck;
%       UPmodel(tt).upids = up;
%       UPmodel(tt).Posterior = Posterior;
%       UPmodel(tt).trials_candidate = trs;
%       UPmodel(tt).Dist = D;
%       UPmodel(tt).target = targets(tr);
%       UPmodel(tt).prediction = C;

      if upcheck && ~isempty(up)
          EP_1block = Epoch_condition(et.Epoch{tr},param_ad);

          %-- update DSP
          EP_update = EP_1block;
          EP_update.nar = EP_update.nar(:,:,up,:);

          %-- update DSP
          param1_ad  = updateDSP(EP_update,C,param1_ad);
          param2_ad  = updateDSP(EP_update,C,param2_ad);

          %-- apply updated DSP to feat

          Repeat = param_ad.repeat;
          param1_ad.repeat = size(EP_update.nar,3);
          param2_ad.repeat = size(EP_update.nar,3);
          param_ad.repeat = size(EP_update.nar,3);
            
          Feat1 = FeatureExt_DSP(EP_update,param1_ad);
          Feat2 = FeatureExt_DSP(EP_update,param2_ad);
          Feat = [Feat1 Feat2];

          [C_new,param_ad] = Classification(Feat,[],param_ad);
 
          label_temp = -ones(size(Feat,1),1);
          label_temp = reshape(label_temp,param1_ad.repeat,param_ad.NumStims);
          label_temp(:,C_new) = 1;
          label = label_temp(:);

          %-- re-calibration
          param_ad.trD.mode = 'training';
          if REMOVE
              

             lbtargets  =  find(param_ad.trD.label == 1);
             lbntargets  =  find(param_ad.trD.label == -1);   

             lbtargetsNew  =  length(find(label == 1));
             lbntargetsNew  =  length(find(label == -1)); 

             IDsNewTarget = lbtargets(lbtargetsNew+1:end);
             IDsNewNTarget = lbntargets(lbntargetsNew+1:end);
             IDsNew = sort([IDsNewTarget;IDsNewNTarget]);

             param_ad.trD.feature = [param_ad.trD.feature(IDsNew,:); Feat];
                param_ad.trD.label = [param_ad.trD.label(IDsNew,:); label];

          else
          param_ad.trD.feature = [param_ad.trD.feature; Feat];
          param_ad.trD.label = [param_ad.trD.label; label];
          end
          
          [~,param_ad] = Classification(param_ad.trD.feature,param_ad.trD.label,param_ad);
          param_ad.repeat = Repeat;
          param1_ad.repeat = Repeat;
          param2_ad.repeat = Repeat;
          %--- 
          %-- re-estimate DSP & SVM
          
%         erptar = EP_update.nar(:,:,:,C);
%         erpnar = EP_update.nar(:,:,:,setdiff(1:param.NumStims,C));
%         X{1} = cat(3,X{1},erptar);
%         X{2} = cat(3,X{2},erpnar(:,:,:));
%         EP_new.tar = X{1};
%         EP_new.nar = X{2};
%         param_ad.trD.mode = 'training';
%         param_ad.DSP = rmfield(param_ad.DSP,'W');
%         [feat_new,label_new, param_ad] = FeatureExt_DSP(EP_new, param_ad);
%         [~,param_ad] = Classification(feat_new,label_new,param_ad);

        %-- 


          fprintf('>> Updated\n')


%           UPmodel(tt).DSP = param_ad.DSP;
%           UPmodel(tt).mdl = param_ad.trD;
%           
%           UPmodel(tt).prediction_new = C_new;
          param_ad.trD.mdl_adapt = param_ad.trD.mdl;
      else
%           UPmodel(tt).prediction_new = NaN;
          fprintf('>> Not updated\n')
      end


        pred_ad(tt) = C;
    
        tt = tt + 1;
    end
    
    pred_ad_con = reshape(pred_ad,Ntr_con,[]);

    Acc_ad = mean(targets_con(:,2:end) == pred_ad_con);
    
if REMOVE
    mkdir([fpath,'\Simulate_',num2str(LAMBDA1),',',num2str(LAMBDA2),'NewRem_Acc'])
    save([fpath,'\Simulate_',num2str(LAMBDA1),',',num2str(LAMBDA2),'NewRem_Acc\',SubName],'Acc_fixed','Acc_ad');
else
    mkdir([fpath,'\Simulate_',num2str(LAMBDA1),',',num2str(LAMBDA2),'New_Acc'])
    save([fpath,'\Simulate_',num2str(LAMBDA1),',',num2str(LAMBDA2),'New_Acc\',SubName],'Acc_fixed','Acc_ad');

end
    % 
%     mkdir([fpath,'\Simulate_',num2str(LAMBDA),'param'])
%     save([fpath,'\Simulate_',num2str(LAMBDA),'param\',SubName],'UPmodel');

end
%%
Acc_fix_all = []; Acc_ad_all = [];
for s =1:Nsub
        SubName = SubNameList{s};
    load([fpath,'\Simulate_',num2str(LAMBDA1),',',num2str(LAMBDA2),'NewRem_Acc\',SubName])
    nsess = size(Acc_fixed,2);
 Acc_fix_all(s,1:nsess) = Acc_fixed;
 Acc_ad_all(s,1:nsess-1) = Acc_ad;
end