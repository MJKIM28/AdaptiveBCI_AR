%%
% compare data selection criteria

% 1) mahalanobis distance vs. 2) 1st 2nd ratio
% (DSP update + Recalibration)

% 1. Accuracy
% 2. Accuracy of included data

clear all; close all

%%


trainIDs = 2:9;
train_order = 2;
prep_use = 'ASRhalf';
lambda = 0.1;
threshold = 0.5;

dataselectiontype = 'ratio';

%%
fpathanalysis = 'D:\[1]EEGBCI\[2]Research\Code&Algorithm\2021_BCI_multitasking\main_2nd\1_Analysis_code';

fpathvars = 'D:\[1]EEGBCI\[1]Data\[11]Multitasking\[2]Main\2nd\2_Behavior\0_Variables';
fpatheeg = 'D:\[1]EEGBCI\[1]Data\[11]Multitasking\[2]Main\2nd\1_EEG';
fpathonline = '\0_online_code\';

files_vars_main = dir([fpathvars,'\*main_var.mat']);
vars_main = {files_vars_main.name};

sublist = cellfun(@(x) x(1:5),vars_main,'UniformOutput',false);

Nsub = length(vars_main);
Ntrial = 25;
Ntrblock = 30;
Ncon  = 3;
Nch = 31;
Nrepeat = 10;
Nrepeat_train = 10;

%%
% prep_use = input('prep type {ASR,ASRcal,ASRhalf,ASR_ICA,ICAfast}:','s');
switch prep_use

    case 'ASR'
        prep_condition = '_ASR';
    case 'ASRcal'
        prep_condition = '_ASRcal';
    case 'ASRhalf'
        prep_condition = '_ASRhalf';
    case 'ASR_ICA'
        prep_condition = '_ASR_ICA';

    case 'ICAfast'
        prep_condition = '_ICAfast';
end

folder_erp = ['ERP',prep_condition];
folder_block = ['Block_data',prep_condition];
folder_param = ['params',prep_condition];
folder_prep = ['Prep',prep_condition];

fpathanalysis_old = fpathanalysis;
fpathparam_m = [fpathanalysis_old,'\Params\'];
fpathtrain = [fpathanalysis_old,'\Train\'];
fpathtest = [fpathanalysis_old,'\Test\'];

fpatherp_train = [fpathtrain,'ERP\',folder_erp];
fpatherp_test = [fpathtest,'ERP\',folder_erp];
fpathparam = [fpathparam_m,folder_param];

cd([pwd,fpathonline])

chanloc = readlocs('D:\[1]EEGBCI\[2]Research\Code&Algorithm\BP_channelLocs.locs');
% train_order = [2 3 1];
%%
global gamma
Acc_ad_both = []; Acc_ad_cl = []; Acc_ad_f = []; Acc_fixed = []; Acc_val = [];
answers_ad_both = []; answers_ad_cl = []; answers_ad_f = []; answers_fixed = [];
Updated = []; Correct = []; Post = []; Dist = [];
%%
for s = 1:Nsub
    subname = sublist{s};
    fprintf('\n\n %s\n',subname);


    e = load([fpatherp_train,'/',subname]);
    te = load([fpatherp_test,'/',subname]);
    load([fpathparam,'/',subname]);

    param = varstosave.param;
    param.trD = rmfield(param.trD,'mdl');
    param.trD.mode = 'training';
    param.DSP = rmfield(param.DSP,'W');
    param.DSP.nf = 10;
    param.winsize = 5;
    param.repeat_init = Nrepeat_train;
    param.repeat = Nrepeat;

    chanlocIn = chanloc;
    chanlocIn(param.badch) = [];
    %%  train

    Ntrained = length(trainIDs);
    EP_train_none = e.varstosave{1};
    EP_train_none.dat = EP_train_none.dat(:,:,1:param.NumStims*param.repeat,trainIDs);
    EP_train_none.lat = EP_train_none.lat(:,1:param.NumStims*param.repeat,trainIDs);
    EP_train_none.target = EP_train_none.target(:,trainIDs);


    [Feature,label,param] = FeatureExtraction_for_analysis(EP_train_none,param);
    param.DSP.init = param.DSP;

    [C,paramL] = Classification(Feature,label,param);
    %     [C,paramL] = Classification_LR(Feature,label,param);

    paramL.trD.feature = Feature;
    paramL.trD.label = label;

    paramNL = paramL;

    paramL.trD.mdl_init = paramL.trD.mdl;
    paramNL.trD.mdl_init = paramNL.trD.mdl;

    %% test (semisupervised online learning)
    %-- new data (one block at a time)
    paramL.repeat = Nrepeat;
    paramNL.repeat = Nrepeat;

    for con = 2
        Ntr = length(te.varstosave{con}.target);

        Thres(s,con) = threshold;
        %%
        %-- validation set (Speak train)
        EP_val = e.varstosave{1};
        EP_val.dat = EP_val.dat(:,:,1:param.NumStims*param.repeat,16:end);
        EP_val.lat = EP_val.lat(:,1:param.NumStims*param.repeat,16:end);
        EP_val.target = EP_val.target(:,16:end);

        EP_val = Epoch_condition(EP_val,paramL);


        for tt = 1:Ntr
            EP_1block = te.varstosave{con};
            EP_1block.dat = EP_1block.dat(:,:,:,tt);
            EP_1block.lat = EP_1block.lat(:,:,tt);
            EP_1block.target = EP_1block.target(tt);

            

            %% get output
            EP_1block = Epoch_condition(EP_1block,paramL);
            
            %% fixed output
            [Feature_te,label_te,paramNL] = FeatureExtraction_newCV(EP_1block,paramNL); %use updated DSP
            [C,paramNL] = Classification(Feature_te,label_te,paramNL); % use updated Classifier

            %% adaptive output
            [Feature_te,label_te,paramL] = FeatureExtraction_newCV(EP_1block,paramL); %use updated DSP
            [~,paramL] = Classification(Feature_te,label_te,paramL); % use updated Classifier
            score = paramL.trD.score;

            %-- initial output (before update)
            output = getoutput(score); %SVM kernel

            %%

            %-- check: use data to update or not

            [up,Posterior,D] = checkupdate(score, [] ,threshold,output,dataselectiontype);

            Updated{tt,con,s} = up;
            Post(:,tt,con,s) = Posterior;
            Dist(:,tt,con,s) = D';

            paramL.DSP.lambda = lambda;

            %% update
            if ~isempty(up)%%&& con ~= 1

                %-- update DSP
                if length(up) > 1
                    EP_update = EP_1block;
                    EP_update.nar = EP_update.nar(:,:,up,:);

                else
                    EP_update = EP_1block;
                end
                paramL  = updateDSP(EP_update,output,paramL);

                %-- apply updated DSP to feat
                paramL.repeat = size(EP_update.nar,3);
                Feat = FeatureExt_DSP(EP_update,paramL);
                [C_new,paramL] = Classification(Feat,[],paramL);

                label_temp = -ones(size(Feat,1),1);
                label_temp = reshape(label_temp,paramL.repeat,param.NumStims);
                label_temp(:,C_new) = 1;
                label = label_temp(:);

                %-- re-calibration
                paramL.trD.feature = [paramL.trD.feature; Feat];
                paramL.trD.label = [paramL.trD.label; label];

                paramL.trD.mode = 'training';
                [~,paramL] = Classification(paramL.trD.feature,paramL.trD.label,paramL);
                %                 paramL = updateClassifier(Feat,output,paramL);
                %                 Nadded = paramL.trD.Nadded;
                %                 Score_update{tt,con,s} = paramL.trD.score_new;

                mdl{tt} = paramL.trD.mdl;
                DSP{tt} = paramL.DSP;

                paramL.repeat = 10;
                paramL.trD.mode = 'testing';

            end


            %% results
            answers_ad_both(tt,con,s) = output;
            answers_fixed(tt,con,s) = C;


        end

        targets(1:Ntr,con,s) = te.varstosave{con}.target';

        Correct(1:Ntr,con,s) =  answers_ad_both(1:Ntr,con,s) == targets(1:Ntr,con,s);

        Acc_ad_both(con,s) = mean(answers_ad_both(1:Ntr,con,s) == targets(1:Ntr,con,s));
        Acc_fixed(con,s) = mean(answers_fixed(1:Ntr,con,s) == targets(1:Ntr,con,s));


    end

end

cd(fpathanalysis)

Correct = logical(Correct);
% savepath = ['D:\[5]defense\4. speech distraction\result_1\'];
% mkdir(savepath);
% savefilename = sprintf(['Accuracy_train%dblock_lambda%.2f_UpAll%d_thres%.1f_trainorder%d%d%d.mat'],Ntrained,paramL.DSP.lambda,UpAll,threshold,train_order);
% save([savepath,savefilename],'Acc_ad_both','Acc_ad_cl','Acc_ad_f','Acc_fixed','Acc_val',...
%     'Correct','Updated','Post','answers_ad_both','answers_ad_cl','answers_ad_f','answers_fixed','targets','Numbers',...
%     'Score_update');


