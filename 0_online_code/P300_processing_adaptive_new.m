function [C, param]   = P300_processing_adaptive_new(sig, trig, param)
global DS
% updated @ 20240311

% try
    if isfield(param.trD, 'mdl')
        param.trD.mode  = 'testing';
        fprintf('Testing start!..\n');
    else
        fprintf('Training start!..\n');
    end

    sig = double(sig*0.04883);
    [sig, param]            = PreProcess(sig,param);

    EP = EpochData(sig,trig,param);
%     [Feature, label, param] = FeatureExtraction(EP, param);
    EP = Epoch_condition(EP,param);
    
    switch param.trD.mode
        case 'training'
            %-- train classifier
            
          [Model_Init, param] = train_model_kfold(param, EP, Nch, Windowlist, Windowlength, Ntime, k);

          param.mdl_init = Model_Init;
%             [C,param]               = Classification(Feature,label, param);
%             param.trD.mdl_init = param.trD.mdl;
%             param.trD.mdl_adapt = param.trD.mdl;
%             param.trD.feature = Feature;
%             param.trD.label = label;
% 
%             param.DSP.init = param.DSP;
        case 'testing'
            
            
             wind = param.Baseline+1:param.Totalepoc;
            EP.nar = preprocess_epochs(EP.nar, wind);
            test_temp = reshape_epochs(EP.nar);

            [rho1, rho2] = compute_rhos(test_temp, param, param.MW.Windowlist, param.MW.Windowlength, param.MW.Ntime);

            [score_1, score_2, pred_1, pred_2, C] = classify_rhos(param.mdl_init.model, rho1, rho2);
            [c1, c2] = get_max_predictions(pred_1, pred_2);
           
            if strcmp(param.trD.ADmode,'fixed')
                C = c2; % output from fixed model
                fprintf('>>> Fixed \n');
                upcheck1 = 0; upcheck2 = 0; 
                D1 = zeros(1,param.NumStims);
                D2 = zeros(1,param.NumStims);
            elseif strcmp(param.trD.ADmode,'adaptive')
                [upcheck1, ~, ~, ~, D1] = checkupdate(score_1, [], [], c1, []);
                [upcheck2, ~, ~, ~, D2] = checkupdate(score_2, [], [], c2, []);
                
                if should_update_model(param.trD.ADmode, param.trD.Fusmode, upcheck1, upcheck2)
                    [param] = update_model(EP, param, param.MW.Ntime, param.trD.Fusmode, c1, c2, upcheck1, upcheck2);
                    fprintf('>>> Block updated\n');
                    
                else
                    fprintf('>>> Block NOT updated\n');
                end
                
            end
            param.updated(1).upcheck = upcheck1;
                param.updated(2).upcheck = upcheck2;
                param.updated(1).D = D1;
                param.updated(2).D = D2;
                param.updated(1).score = score_1;
                param.updated(2).score = score_2;
 
%             if strcmp(param.trD.ADmode,'fixed')
%                 param.trD.mdl = param.trD.mdl_init;
%             [C,param] = Classification(Feature,label,param); % use updated Classifier
%                      fprintf('>> Use fixed\n') 
%       
% 
%             elseif strcmp(param.trD.ADmode,'adaptive')
%                 param.trD.mdl = param.trD.mdl_adapt;
%                  %-- initial output
%             [C,param] = Classification(Feature,label,param); % use updated Classifier
%             score = param.trD.score;
%             output = getoutput(score); %SVM kernel
           
%                 %-- check update
%                 [upcheck,up,Posterior,trs,D] = checkupdate(score, [] ,param.trD.threshold,output,[]);
%                 if upcheck
%                     param.update{param.Numtrial}.updated = up;
%                 else
%                     param.update{param.Numtrial}.updated = 0;
%                 end
%                 param.update{param.Numtrial}.Posterior = Posterior;
%                 param.update{param.Numtrial}.Dist = D;
% 
%                 %-- update
%                 %-- DSP
%                 if upcheck
%                     EP_1block = Epoch_condition(EP,param);
% 
%                     %-- update DSP
%                     EP_update = EP_1block;
%                     EP_update.nar = EP_update.nar(:,:,up,:);
% 
%                     param  = updateDSP(EP_update,output,param);
% 
%                     %-- apply updated DSP to feat
% 
%                     Repeat = param.repeat;
%                     param.repeat = size(EP_update.nar,3);
%                     Feat = FeatureExt_DSP(EP_update,param);
%                     [C_new,param] = Classification(Feat,[],param);
% 
%                     label_temp = -ones(size(Feat,1),1);
%                     label_temp = reshape(label_temp,param.repeat,param.NumStims);
%                     label_temp(:,C_new) = 1;
%                     label = label_temp(:);
% 
%                     %-- re-calibration
%                     param.trD.feature = [param.trD.feature; Feat];
%                     param.trD.label = [param.trD.label; label];
% 
%                     param.trD.mode = 'training';
%                     [~,param] = Classification(param.trD.feature,param.trD.label,param);
% 
%                     param.repeat = Repeat;
% 
%                     fprintf('>> Updated\n')
% 
% 
%                     param.update{param.Numtrial}.DSP = param.DSP;
%                     param.update{param.Numtrial}.mdl = param.trD;
%                     param.trD.mdl_adapt = param.trD.mdl;
%                 else
%                     fprintf('>> Not updated\n')
%                 end          
%             end

            % FOR Dynamic Stopping
            % If iteration end but block not end,
            % decide to stop or not
            if (~param.switch_on) && isfield(param,'switch_on_iter')
                dat = reshape(param.decoder.data,param.NumStims,size(param.decoder.data,1)/param.NumStims);

                DS.data = dat;
                DS = getpval(DS);
                [param.DS.stop,param.DS.class] = decidestopping(DS);
            end
    end
    % clearvars -except C param
% catch
%     keyboard
% end
end

function [Model_Init, param_temp] = train_model_kfold(EP,param)

NFmax = fix(param.NumCh/ 2);
Acc_val = zeros(NFmax, param.MW.Ntime, param.MW.kfold);


param.trD.mode = 'training';


Nblock =size(EP.target, 2);
Nbperfold = round(Nblock/k);
rng(1)
folds = crossvalind('Kfold', Nblock, k);

Nsample = Nbperfold*param.repeat*param.NumStims;
all_yy = zeros(Nsample,param.MW.Ntime,NFmax,k);
all_lbl = zeros(Nsample,param.MW.Ntime,NFmax,k);
wind = param.Baseline + 1:param.Totalepoc;
fprintf('k fold')
for fold = 1:k
    train_idx = folds ~= fold;
    val_idx = folds == fold;


    EP_train_tr_fold = select_epochs(EP, train_idx,param.repeat,wind);
    EP_train_val_fold = select_epochs(EP, val_idx,param.repeat,wind);

    [val_temp_fold, Label_val] = prepare_validation_data(EP_train_val_fold);

    [Acc,yy, lbl] = validate_NF(param, EP_train_tr_fold, val_temp_fold, Label_val, NFmax, param.MW.Windowlist, param.MW.Windowlength, param.MW.Ntime);
    Acc_val(:,:,fold) = Acc;
    all_yy(1:size(yy,1),:,:,fold) = yy;
    all_lbl(1:size(lbl,1),:,:,fold) =lbl;
    fprintf('.')
end
fprintf('\n')
[~, Nfs] = max(mean(Acc_val, 3));  % Mean accuracy across folds

yy_select = zeros(Nsample,k,Ntime);
lbl_select = zeros(Nsample,k,Ntime);
for t = 1:Ntime
    yy_select(:,:,t) = squeeze(all_yy(:,t,Nfs(t),:));
    lbl_select(:,:,t) = squeeze(all_lbl(:,t,Nfs(t),:));
end
yy_select = reshape(yy_select,[],Ntime);
lbl_select = reshape(lbl_select,[],Ntime);

% Train logistic regression on all folds' data
model = fitglm(yy_select, lbl_select(:,1), 'linear', 'Distribution', 'binomial');


% Initialize final model with all data
EP_train_all = select_epochs(EP, [1:Nblock]',param.repeat,wind);

[Model_Init, param_temp] = initialize_model(param, EP_train_all, Nfs, param.MW.Windowlist, param.MW.Windowlength, param.MW.Ntime, model);
end


function EP_selected = select_epochs(EP,idx,Nrepeat,wind)
EP_selected.target = EP.target(idx);
idx_l = find(idx);
tridall = bsxfun(@plus, (idx_l-1) * Nrepeat, [1:Nrepeat])';
tridall2 = bsxfun(@plus, (idx_l-1) * Nrepeat * 3, [1:Nrepeat * 3])';
EP_selected.tar = EP.tar(wind, :, tridall(:));
EP_selected.nar = EP.nar(wind, :, tridall2(:));


EP_selected.tar = EP_selected.tar - mean(EP_selected.tar, 1, 'omitnan');
EP_selected.nar = EP_selected.nar - mean(EP_selected.nar, 1, 'omitnan');
end

function [val_temp, Label_val, Nval] = prepare_validation_data(EP_train_val)
val_temp = cat(3, EP_train_val.tar, EP_train_val.nar);
val_temp = permute(val_temp, [2, 1, 3]);
Nval = size(val_temp, 3);
Label_val = [ones(size(EP_train_val.tar, 3), 1); zeros(size(EP_train_val.nar, 3), 1)];
end


function [Acc_val,yy, lbl] = validate_NF(param, EP_train_tr, val_temp, Label_val, NFmax, Windowlist, Windowlength, Ntime)
Acc_val = zeros(NFmax,Ntime);
yy = zeros(length(Label_val),Ntime,NFmax);
lbl = repmat(Label_val,1,Ntime,NFmax);
for Nf = 1:NFmax
    param.DSP.nf = Nf;

    for ti = 1:Ntime
        [DSP1, LDA1] = get_window_models(param, EP_train_tr, Windowlist(ti), Windowlength, Nf);

        [X_,yy_temp] = extract_features_n_outputs(DSP1, val_temp, LDA1, DSP1.window);
        yy(:,ti,Nf) = yy_temp;
    end
    Acc_val(Nf,:) = sum((yy(:,:,Nf) > 0.5) == Label_val) / length(Label_val);

end
end

function [DSP1, LDA1] = get_window_models(param, EP_train_tr, window_start, Windowlength, Nf)
windowind = window_start:5:(window_start + Windowlength - 1);
param.DSP.window = windowind;
param.DSP.nf = Nf;

param.trD.mode = 'training';

[~, ~, param_temp] = FeatureExt_DSP(EP_train_tr, param);
[Template_Tar, Template_Nar] = filter_templates(param_temp.DSP, EP_train_tr, windowind, Nf);
X = [reshape(Template_Tar, [], size(Template_Tar, 3))'; reshape(Template_Nar, [], size(Template_Nar, 3))'];
Y = [ones(size(Template_Tar, 3), 1); zeros(size(Template_Nar, 3), 1)];
DSP1 = param_temp.DSP;
LDA1 = getLDA(X, Y);

end

function [Template_Tar, Template_Nar] = filter_templates(DSP1, EP_train_tr, windowind, Nf)
Ntr = size(EP_train_tr.tar,3);
Ntr_nar = size(EP_train_tr.nar,3);
Nt = length(windowind);
Template_Tar = zeros(Nf,Nt,Ntr);
Template_Nar = zeros(Nf,Nt,Ntr_nar);
for tr = 1:Ntr
    Template_Tar(:, :, tr) = DSP1.W(:, 1:Nf)' * EP_train_tr.tar(windowind, :, tr)';
end
for tr = 1:Ntr_nar
    Template_Nar(:, :, tr) = DSP1.W(:, 1:Nf)' * EP_train_tr.nar(windowind, :, tr)';
end
end


function [X_,yy] = extract_features_n_outputs(DSP1, val_temp, mdl_mj, windowind)
Ntr = size(val_temp,3);
Nt = length(windowind);
X_ = zeros(DSP1.nf,Nt,Ntr);
for tr = 1:size(val_temp, 3)
    X_(:, :, tr) = DSP1.W(:, 1:DSP1.nf)' * val_temp(:, windowind, tr);
end
[sc_2] = pred_mj(mdl_mj, reshape(X_, [], size(X_, 3))');
yy = sc_2(:, end);
end


function [Model_Init, param_new] = initialize_model(param, EP_train_tr,  Nfs, Windowlist, Windowlength, Ntime, model)
% train from scratch

param_new = param;
if isfield(param_new,'DSP')
    param_new = rmfield(param_new,'DSP');
end
if isfield(param_new.trD,'mdl')
    param_new.trD = rmfield(param_new.trD,'mdl');
end
for ti = 1:Ntime
    [DSP1, LDA1] = get_window_models(param, EP_train_tr, Windowlist(ti), Windowlength, Nfs(ti));
    LDA1.lambda = 0.03;
    LDA1.CosSim = NaN;
    DSP1.lambda = 0.03;
    DSP1.CosSim = NaN(param.NumCh,1);

    param_new.DSP(ti) = DSP1;
    param_new.trD.mdl(ti) = LDA1;
end
% param_new.trD.lastdecisionmdl = model;

Model_Init.DSP = param_new.DSP;
Model_Init.LDA = param_new.trD.mdl;
Model_Init.model = model;

end



function nar = preprocess_epochs(nar, wind)
nar = nar(wind, :, :, :);
nar = nar - mean(nar, 1, 'omitnan');
end

function reshaped = reshape_epochs(nar)
reshaped = permute(nar, [2, 1, 3, 4]);
reshaped = reshape(reshaped, size(reshaped, 1), size(reshaped, 2), []);
end

function [rho1, rho2] = compute_rhos(test_temp, param,  Windowlist_new, Windowlength, Ntime_new)
rho1 = NaN(size(test_temp, 3), Ntime_new);
rho2 = NaN(size(test_temp, 3), Ntime_new);

for ti = 1:Ntime_new
    windowind = Windowlist_new(ti):5:Windowlist_new(ti) + Windowlength - 1;

    W1 = param.DSP(ti).W(:, 1 : param.DSP(ti).nf);
    W2 = param.mdl_init.DSP(ti).W(:, 1 : param.mdl_init.DSP(ti).nf);

    for tt = 1:size(test_temp, 3)
        Y_1 = W1' * test_temp(:, windowind, tt);
        Y_2 = W2' * test_temp(:, windowind, tt);

        [sc_11] = pred_mj(param.trD.mdl(ti), Y_1(:)');
        [sc_112] = pred_mj(param.mdl_init.LDA(ti), Y_2(:)');

        rho1(tt, ti) = sc_11(:, end);
        rho2(tt, ti) = sc_112(:, end);
    end
end
end

function [score_1, score_2, pred_1, pred_2, predicted] = classify_rhos(model, rho1, rho2)
sc2 = predict(model, rho1);
score_1 = reshape(sum(sc2, 2), 10, 4);
pred_1 = sum(score_1);

sc22 = predict(model, rho2);
score_2 = reshape(sum(sc22, 2), 10, 4);
pred_2 = sum(score_2);

[~, predicted] = max(normalize(pred_1) + normalize(pred_2));
end

function [c1, c2] = get_max_predictions(pred_1, pred_2)
[~, c1] = max(pred_1);
[~, c2] = max(pred_2);
end


function decision = should_update_model(adaptmode, fusionmode, upcheck1, upcheck2)
if strcmp(adaptmode, 'adaptive')
    if strcmp(fusionmode, 'fusion')
        decision = upcheck1 || upcheck2;
    else
        decision = upcheck1;
    end
else
    decision = false;
end
end

function [param] = update_model(EP_update, param, Ntime_new, fusionmode, c1, c2, upcheck1, upcheck2)
if strcmp(fusionmode, 'fusion')
    if upcheck1
        uplabel = c1;
    elseif upcheck2
        uplabel = c2;
    end
else
    uplabel = c1;
end


% yynew = [];
for ti = 1:Ntime_new
    param_temp = param;
    param_temp.DSP  = param.DSP(ti);
    param_temp = updateDSP(EP_update, uplabel, param_temp);
    param.DSP(ti) = param_temp.DSP;

    Feat = FeatureExt_DSP(EP_update, param_temp);
    LDA = updateSingleClassfier(Feat, param.trD.mdl(ti), []);
    param.trD.mdl(ti) = LDA;
%     [sc_2] = pred_mj(param.trD.mdl(ti), Feat);
%     yynew(:, ti) = sc_2(:, end);
end
end