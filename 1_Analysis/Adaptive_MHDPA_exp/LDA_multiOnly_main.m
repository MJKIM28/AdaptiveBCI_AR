
function LDA_multiOnly_main(adaptmode,fusionmode)
fpath = 'E:\1_연구\분석코드\2024_adaptive_BCI\Data\Pilot';

fpathanalysis = cd;



%%

SubNames = {'Subtest12','Subtest13','Subtest14','Subtest15','Subtest17',...
    'Subtest18','Subtest19','Subtest20','Subtest21','Subtest22'};
Nsub = length(SubNames);

fpathonline ='E:\1_연구\분석코드\2024_adaptive_BCI\AdaptiveBCI_AR\0_online_code';


CondName = {'Pre','Main 1','Main 2','Main 3','Main 4','Main 5','Post'};
Ncon = length(CondName);

Ntrial = 15;
Nrepeat = 10;

ChNames = {'Fp1'    'Fpz'    'Fp2'    'F7'    'F3'    'Fz'    'F4'    'F8' ...
    'FT9'    'FC5'    'FC1'    'FC2'    'FC6'    'FT10'    'T7'    'C3'...
    'Cz'    'C4'    'T8'    'CP5'    'CP1'    'CP2'    'CP6'    'P7'...
    'P3'    'Pz'    'P4'    'P8'    'O1'    'Oz'    'O2'};
Nch = length(ChNames);

Fs = 500;
Time = -0.2:1/Fs:0.6-1/Fs;
Time2 = 0:1/Fs:0.6-1/Fs;

Time_d = Time(1:10:end);



Windowlength = fix(0.2*Fs);
Windowlist = fix(1:0.05*Fs:0.6*Fs-Windowlength);
Ntime = length(Windowlist);

chanloc = readlocs('D:\[1]EEGBCI\[2]Research\Code&Algorithm\BP_channelLocs.locs');

load('E:\1_연구\분석코드\2024_adaptive_BCI\AdaptiveBCI_AR\1_Analysis\TrialUsed2.mat')
load('E:\1_연구\분석코드\2024_adaptive_BCI\AdaptiveBCI_AR\1_Analysis\TrialUsed.mat')


trial_use1 = Trials;
trial_use1 = trial_use1(9:end);
trial_use2 = Trials2;
trial_use = [trial_use1;trial_use2];


scores1 = NaN(40,Ntime,Ntrial,Ncon,Nsub);
scores2 = NaN(40,Ntime,Ntrial,Ncon,Nsub);
outputs1 = NaN(4,Ntrial,Ncon,Nsub);
outputs2 = NaN(4,Ntrial,Ncon,Nsub);

outputs = NaN(Ntrial,Ncon,Nsub);
answers = NaN(Ntrial,Ncon,Nsub);
DSPs = [];
        cd(fpathonline)

for s = 1:Nsub
    
if s <5
    Ntrial = 15;
else
    Ntrial = 14;
end

    SubName = SubNames{s};
    fprintf('%s\n',SubName);

    e = load([fpath,'\Train\Epoch\',SubName]);
    et = load([fpath,'\Test\Epoch\',SubName]);
    p = load([fpath,'\Dat_',SubName,'\param.mat']);
    param = p.param;

    chanlocIn = chanloc;
    chanlocIn(param.badch) = [];

    badch = param.badch;
    chlist = setdiff(1:Nch,badch);

    param.DSP = rmfield(param.DSP,'W');
    %%



    param.trD.mode = 'training';
    param.repeat = Nrepeat;
    EP = e.Epoch;
   
    %-- only control 15 blocks
    trid = [1:10];
    valid = [11:15];
    EP_train_tr.target = EP.target(trid);
    tridall = bsxfun(@plus,(trid-1)*Nrepeat,[1:Nrepeat]');
    tridall2 = bsxfun(@plus,(trid-1)*Nrepeat*3,[1:Nrepeat*3]');
    EP_train_tr.tar = EP.tar(0.2*Fs+1:0.8*Fs,:,tridall(:));
    EP_train_tr.nar = EP.nar(0.2*Fs+1:0.8*Fs,:,tridall2(:));

    EP_train_val.target = EP.target(valid);
    validall = bsxfun(@plus,(valid-1)*Nrepeat,[1:Nrepeat]');
    validall2 = bsxfun(@plus,(valid-1)*Nrepeat*3,[1:Nrepeat*3]');
    EP_train_val.tar = EP.tar(0.2*Fs+1:0.8*Fs,:,validall);
    EP_train_val.nar = EP.nar(0.2*Fs+1:0.8*Fs,:,validall2);


    EP_train_tr.tar = EP_train_tr.tar - mean(EP_train_tr.tar,1,'omitnan');
    EP_train_tr.nar = EP_train_tr.nar - mean(EP_train_tr.nar,1,'omitnan');
    EP_train_val.tar = EP_train_val.tar - mean(EP_train_val.tar,1,'omitnan');
    EP_train_val.nar = EP_train_val.nar - mean(EP_train_val.nar,1,'omitnan');

    EP_all.target = [EP_train_tr.target EP_train_val.target];
    EP_all.tar = cat(3,EP_train_tr.tar(1:5:end,:,:),EP_train_val.tar(1:5:end,:,:));
    EP_all.nar = cat(3,EP_train_tr.nar(1:5:end,:,:),EP_train_val.nar(1:5:end,:,:));

    val_temp = cat(3,EP_train_val.tar,EP_train_val.nar);
    val_temp = permute(val_temp,[2,1,3]);
    Nval = size(val_temp,3);
    Label_val = [ones(size(EP_train_val.tar,3),1);zeros(size(EP_train_val.nar,3),1)];


% validation - NF
    for Nf = 1:fix(length(chlist)/2)
        param.DSP.nf = Nf;
        yy = []; mdl = []; param_temp = [];
        for ti = 1:Ntime
            windowind = Windowlist(ti):5:Windowlist(ti)+Windowlength-1;
            param.DSP.window = windowind;

            EP_input = EP_train_tr;
            [~,~,param_temp{ti}] = FeatureExt_DSP(EP_input,param);


            %-- filtered template
            Template_Tar = []; Template_Nar = [];
            for tr = 1:size(EP_train_tr.tar,3)
                Template_Tar(:,:,tr) = param_temp{ti}.DSP.W(:,1: param.DSP.nf)'*EP_train_tr.tar(windowind,:,tr)';
            end
            for tr = 1:size(EP_train_tr.nar,3)
                Template_Nar(:,:,tr) = param_temp{ti}.DSP.W(:,1: param.DSP.nf)'*EP_train_tr.nar(windowind,:,tr)';
            end

            X = [reshape(Template_Tar,[],size(Template_Tar,3))';reshape(Template_Nar,[],size(Template_Nar,3))'];
            Y = [ones(size(Template_Tar,3),1);zeros(size(Template_Nar,3),1)];

            mdl_mj{ti} = getLDA(X,Y);
            % validation
            X_ = [];
            for tr = 1:Nval
                X_(:,:,tr) = param_temp{ti}.DSP.W(:,1: param.DSP.nf)'*val_temp(:,windowind,tr);
            end

            [sc_2] = pred_mj(mdl_mj{ti},reshape(X_,[],size(X_,3))');

            yy(:,ti) = sc_2(:,end);

        end
        Acc_val(Nf,:) = sum((yy>0.5) == Label_val)/length(Label_val);

    end

     [~,id] = max(Acc_val);
    Nfs = id;


    %--- TRAIN   
    yy = []; mdl = []; param_temp = [];
    for ti = 1:Ntime
        windowind = Windowlist(ti):5:Windowlist(ti)+Windowlength-1;
        param.DSP.window = windowind;
        param.DSP.nf = Nfs(ti);
        
        EP_input = EP_train_tr;
        [~,~,param_temp{ti}] = FeatureExt_DSP(EP_input,param);
        
        %-- filtered template
        Template_Tar = []; Template_Nar = [];
        for tr = 1:size(EP_train_tr.tar,3)
            Template_Tar(:,:,tr) = param_temp{ti}.DSP.W(:,1: param.DSP.nf)'*EP_train_tr.tar(windowind,:,tr)';
        end
        for tr = 1:size(EP_train_tr.nar,3)
            Template_Nar(:,:,tr) = param_temp{ti}.DSP.W(:,1: param.DSP.nf)'*EP_train_tr.nar(windowind,:,tr)';
        end
        
        X = [reshape(Template_Tar,[],size(Template_Tar,3))';reshape(Template_Nar,[],size(Template_Nar,3))'];
        Y = [ones(size(Template_Tar,3),1);zeros(size(Template_Nar,3),1)];
        
        mdl_mj{ti} = getLDA(X,Y);
      
        % validation
        X_ = [];
        for tr = 1:Nval
            X_(:,:,tr) = param_temp{ti}.DSP.W(:,1: param.DSP.nf)'*val_temp(:,windowind,tr);
        end
        [sc_2] = pred_mj(mdl_mj{ti},reshape(X_,[],size(X_,3))');
        
        yy(:,ti) = sc_2(:,end);
        
        param_temp{ti}.DSP.lambda = NaN;
        param_temp{ti}.trD.mode = 'testing';
        mdl_mj{ti}.lambda = NaN;
        
        param_temp{ti}.DSP.CosSim = NaN(length(chlist),1);
        mdl_mj{ti}.CosSim = NaN;
    end



    % % 로지스틱 회귀 모델 초기 훈련
    X_init = yy; Y_init = Label_val;
    model = fitglm(X_init, Y_init, 'linear', 'Distribution', 'binomial');

    % 계수의 절대값을 사용하여 중요도 평가

    inmodel = 1:size(X_init,2);

    Windowlist_new = Windowlist(inmodel);
    Ntime_new = length(Windowlist_new);


Model_Init.DSP = param_temp;
Model_Init.LDA = mdl_mj;
    BadChs{s} = badch;
%     Temporal{s}.estimate = model.Coefficients.Estimate;
%     Temporal{s}.pval = model.Coefficients.pValue;
%     DSPtotal{s} = param_all.DSP.W;
   LastDecision{s} = model;

    Xnew = X_init;
    Ynew = Y_init;
   %%
    for c = 1:Ncon
        %--- test block
       
        trials = find(trial_use{s} == Ntrial*(c-1)+1):find(trial_use{s} == Ntrial*c);
        %%
if ~isempty(trials)
     t = 1;
        for tr = trials
            EP_T = et.Epoch{tr};

            EP_Te.dat = EP_T.dat(:,:,1:4*Nrepeat);
            EP_Te.lat = EP_T.lat(:,1:4*Nrepeat);


            param.trD.mode = 'testing';
            EP_test = Epoch_condition(EP_Te,param);

            %--- initial classification
            EP_test.nar = EP_test.nar(0.2*Fs+1:0.8*Fs,:,:,:);
            EP_test.nar = EP_test.nar - mean(EP_test.nar,1,'omitnan');

            test_temp = permute(EP_test.nar,[2,1,3,4]);
            test_temp = reshape(test_temp,size(test_temp,1),size(test_temp,2),[]);

            rho1 = []; Feat_te = []; rho2 = [];
            for tt = 1:size(test_temp,3)
                for ti = 1:Ntime_new

                    windowind = Windowlist_new(ti):5:Windowlist_new(ti)+Windowlength-1;

                    Y_ = param_temp{ti}.DSP.W(:,1: param_temp{ti}.DSP.nf)'*test_temp(:,windowind,tt);
                    [sc_11] = pred_mj(mdl_mj{ti},Y_(:)');
                    rho1(tt,ti) = sc_11(:,end);


                       % fixed
                        Y_ = Model_Init.DSP{ti}.DSP.W(:,1: param_temp{ti}.DSP.nf)'*test_temp(:,windowind,tt);
                        [sc_112] = pred_mj(Model_Init.LDA{ti},Y_(:)');
                        rho2(tt,ti) = sc_112(:,end);

                end

            end
            [sc2] = predict(model,rho1);
            score_1 = reshape(sum(sc2,2),10,4);
            pred_1 =  sum(score_1);

            [sc22] = predict(model,rho2);
            score_2 = reshape(sum(sc22,2),10,4);
            pred_2 =  sum(score_2);


            [~,predicted] = max(normalize(pred_1)+normalize(pred_2));

            if c > 1

            [~,c1] = max(pred_1);
            [~,c2] = max(pred_2);

            %--- data selection
            [upcheck1,up1,Posterior1,trs1,D1] = checkupdate(score_1, [] ,[],c1,[]);
            [upcheck2,up2,Posterior2,trs2,D2] = checkupdate(score_2, [] ,[],c2,[]);

%             w1 = 1; w2 = 1;
%            if strcmp(fusionmode,'fusion')
%                if upcheck1 
%                    w1 = 1.5;
%                end
%                if upcheck2
%                    w2 = 1.5;
%                end
%                [~,predicted] = max(normalize(pred_1)*w1+normalize(pred_2)*w2);
%                
% 
%            end


            %--- update model
            if strcmp(adaptmode,'adaptive')

                if strcmp(fusionmode,'fusion')
                    decision = upcheck1 || upcheck2;
                else
                    decision = upcheck1;
                end
                if decision
                    if strcmp(fusionmode,'fusion')
                        if upcheck1
                            uplabel = c1;

                        elseif upcheck2
                            uplabel = c2;

                        end
                    else
                        uplabel = c1;
                    end
                fprintf('Block %d updated\n',tr)

                %-- update DSP
                EP_update = EP_test;
                Feat = []; yynew = [];
                for ti = 1:Ntime_new

                    param_temp{ti}  = updateDSP(EP_update,uplabel,param_temp{ti});
                    Feat{ti} = FeatureExt_DSP(EP_update,param_temp{ti});

                    mdl_mj{ti} = updateSingleClassfier(Feat{ti},mdl_mj{ti},[]);


                    [sc_2] = pred_mj(mdl_mj{ti},Feat{ti});

                    yynew(:,ti) = sc_2(:,end);

                end

                if size(Xnew,1) < Ntime*50
                label_temp = zeros(size(Feat{1},1),1);
                label_temp = reshape(label_temp,Nrepeat,4);
                label_temp(:,uplabel) = 1;
                label = label_temp(:);
                Xnew = [Xnew; yynew]; Ynew = [Ynew; label]; 

                model = fitglm(Xnew, Ynew, 'linear', 'Distribution', 'binomial');
                end


            else
                fprintf('Block %d NOT updated\n',tr);
                end
            end

            dsptemp = [];
            for ti = 1:Ntime_new
                dsptemp(:,:,ti) = param_temp{ti}.DSP.W;
                lambdas(ti) = param_temp{ti}.DSP.lambda;
            end
            DSPs{t,c,s} = dsptemp;
            Lambdas{t,c,s} = lambdas;
            LDA{t,c,s} = mdl_mj;
            UpDistance1(:,t,c,s) = D1;
            UpDistance2(:,t,c,s) = D2;
            Updated1(t,c,s) = upcheck1;
            Updated2(t,c,s) = upcheck2;

            end


            scores1(:,:,t,c,s) = rho1;
            scores2(:,:,t,c,s) = rho2;
            outputs1(:,t,c,s) = pred_1;
            outputs2(:,t,c,s) = pred_2;

            outputs(t,c,s) = predicted;
            answers(t,c,s) = et.Epoch{tr}.target;

  
 t = t+1;
        end
       

end
    end
end
cd(fpathanalysis)
mkdir('adaptiveMWHDPA_LDA_MultiOnly_NoThre_VarNf_VarUC_updateLRmax50')
save(['adaptiveMWHDPA_LDA_MultiOnly_NoThre_VarNf_VarUC_updateLRmax50\',adaptmode,'_',fusionmode,'.mat'],...
     'scores1','scores2','outputs1','outputs2','outputs','answers',...
     'DSPs','BadChs','LDA','Updated1','Updated2','UpDistance1','UpDistance2','Lambdas','LastDecision')
