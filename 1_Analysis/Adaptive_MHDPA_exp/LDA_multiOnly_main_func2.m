function LDA_multiOnly_main_func2(adaptmode, fusionmode, k)
if nargin < 3
    k = 5; % Default number of folds if not provided
end

% Set file paths
fpath = 'E:\1_연구\분석코드\2024_adaptive_BCI\Data\Pilot';
fpathonline = 'E:\1_연구\분석코드\2024_adaptive_BCI\AdaptiveBCI_AR\0_online_code';
fpathanalysis = cd;

% Initialize variables
SubNames = {'Subtest12', 'Subtest13', 'Subtest14', 'Subtest15', 'Subtest17', ...
    'Subtest18', 'Subtest19', 'Subtest20', 'Subtest21', 'Subtest22'};
Nsub = length(SubNames);
CondName = {'Pre', 'Main 1', 'Main 2', 'Main 3', 'Main 4', 'Main 5', 'Post'};
Ncon = length(CondName);
Nrepeat = 10;
Fs = 500;
Windowlength = fix(0.2 * Fs);
Windowlist = fix(1: 0.05 * Fs : 0.6 * Fs - Windowlength);
Ntime = length(Windowlist);
chanloc = readlocs('D:\[1]EEGBCI\[2]Research\Code&Algorithm\BP_channelLocs.locs');
trial_use = load_trials();

% Initialize storage variables
Ntrial = 15;
[scores1, scores2, outputs1, outputs2, outputs, answers, DSPs, BadChs, Lambdas, LastDecision, UpDistance1, UpDistance2, Updated1, Updated2, LDA] = initialize_storage(Ntime, Ntrial, Ncon, Nsub);
mdl_init = cell(Nsub,1);
for s = 1:Nsub

    if s <5
    Ntrial = 15;
else
    Ntrial = 14;
    end

    SubName = SubNames{s};
    fprintf('%s\n', SubName);

    [e, et, param, chlist] = load_subject_data(fpath, SubName, chanloc,Nrepeat);
    %         [param, EP_train_tr, EP_train_val, val_temp, Label_val, Nval] = prepare_training_data(e, param, Nrepeat, Fs);

    % Training Part with k-fold cross-validation
    Nch = length(chlist);
    [Model_Init, param_temp] = train_model_kfold(param, e, Nch, Windowlist, Windowlength, Ntime, k);
    mdl_init{s} = Model_Init;
    % Testing Part
    [scores1, scores2, outputs1, outputs2, outputs, answers, DSPs, Lambdas, LDA, UpDistance1, UpDistance2, Updated1, Updated2] = ...
        test_blocks(adaptmode, fusionmode, Ncon, Ntrial, trial_use, s, Nrepeat, et, param_temp, Model_Init, Ntime, Windowlist, Windowlength, ...
        scores1, scores2, outputs1, outputs2, outputs, answers, DSPs, Lambdas, LDA, UpDistance1, UpDistance2, Updated1, Updated2);
end

save_results(fpathanalysis, adaptmode, fusionmode, scores1, scores2, outputs1, outputs2, outputs, answers, DSPs, BadChs, LDA, Updated1, Updated2, UpDistance1, UpDistance2, Lambdas, mdl_init);
end



function trial_use = load_trials()
load('E:\1_연구\분석코드\2024_adaptive_BCI\AdaptiveBCI_AR\1_Analysis\TrialUsed2.mat');
load('E:\1_연구\분석코드\2024_adaptive_BCI\AdaptiveBCI_AR\1_Analysis\TrialUsed.mat');
trial_use1 = Trials;
trial_use1 = trial_use1(9:end);
trial_use2 = Trials2;
trial_use = [trial_use1; trial_use2];
end

function [scores1, scores2, outputs1, outputs2, outputs, answers, DSPs, BadChs, Lambdas, LastDecision, UpDistance1, UpDistance2, Updated1, Updated2, LDA] = initialize_storage(Ntime, Ntrial, Ncon, Nsub)
scores1 = NaN(40, Ntime, Ntrial, Ncon, Nsub);
scores2 = NaN(40, Ntime, Ntrial, Ncon, Nsub);
outputs1 = NaN(4, Ntrial, Ncon, Nsub);
outputs2 = NaN(4, Ntrial, Ncon, Nsub);
outputs = NaN(Ntrial, Ncon, Nsub);
answers = NaN(Ntrial, Ncon, Nsub);
DSPs = cell(Ntrial,Ncon,Nsub);
BadChs = cell(Nsub, 1);
Lambdas = NaN(Ntime,Ntrial, Ncon, Nsub);
LDA = cell(Ntrial, Ncon, Nsub);
Updated1 = NaN(Ntrial, Ncon, Nsub);
Updated2 = NaN(Ntrial, Ncon, Nsub);
UpDistance1 = NaN(4,Ntrial, Ncon, Nsub);
UpDistance2 = NaN(4,Ntrial, Ncon, Nsub);
LastDecision = cell(Nsub, 1);
end

function [e, et, param, chlist] = load_subject_data(fpath, SubName, chanloc, Nrepeat)
e = load([fpath, '\Train\Epoch\', SubName]);
et = load([fpath, '\Test\Epoch\', SubName]);
p = load([fpath, '\Dat_', SubName, '\param.mat']);
param = p.param;
param.DSP = rmfield(param.DSP, 'W');
param.chanlocIn = chanloc;
param.chanlocIn(param.badch) = [];
chlist = setdiff(1:length(chanloc), param.badch);
param.repeat = Nrepeat;
end


function [Model_Init, param_temp] = train_model_kfold(param, e, Nch, Windowlist, Windowlength, Ntime, k)

NFmax = fix(Nch / 2);
Acc_val = zeros(NFmax, Ntime, k);


param.trD.mode = 'training';
EP = e.Epoch;


Nblock =size(EP.target, 2);
Nbperfold = fix(Nblock/k);
rng(1)
folds = crossvalind('Kfold', Nblock, k);

Nsample = Nbperfold*param.repeat*param.NumStims;
all_yy = zeros(Nsample,Ntime,NFmax,k);
all_lbl = zeros(Nsample,Ntime,NFmax,k);
wind = param.Baseline + 1:param.Totalepoc;
fprintf('k fold')
for fold = 1:k
    train_idx = folds ~= fold;
    val_idx = folds == fold;


    EP_train_tr_fold = select_epochs(EP, train_idx,param.repeat,wind);
    EP_train_val_fold = select_epochs(EP, val_idx,param.repeat,wind);

    [val_temp_fold, Label_val] = prepare_validation_data(EP_train_val_fold);

    [Acc,yy, lbl] = validate_NF(param, EP_train_tr_fold, val_temp_fold, Label_val, NFmax, Windowlist, Windowlength, Ntime);
    Acc_val(:,:,fold) = Acc;
    all_yy(:,:,:,fold) = yy;
    all_lbl(:,:,:,fold) =lbl;
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

[Model_Init, param_temp] = initialize_model(param, EP_train_all, Nfs, Windowlist, Windowlength, Ntime, model);
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
    LDA1.CosSim = [];
    DSP1.lambda = 0.03;
    DSP1.CosSim = [];

    param_new.DSP(ti) = DSP1;
    param_new.trD.mdl(ti) = LDA1;
end
param_new.trD.lastdecisionmdl = model;

Model_Init.DSP = param_new.DSP;
Model_Init.LDA = param_new.trD.mdl;
Model_Init.model = model;

end


function [scores1, scores2, outputs1, outputs2, outputs, answers, DSPs, Lambdas, LDA, UpDistance1, UpDistance2, Updated1, Updated2] = test_blocks(adaptmode, fusionmode, Ncon, Ntrial, trial_use, s, Nrepeat, et, param_temp, Model_Init, Ntime_new, Windowlist_new, Windowlength,scores1, scores2, outputs1, outputs2, outputs, answers, DSPs, Lambdas, LDA, UpDistance1, UpDistance2, Updated1, Updated2)
for c = 1:Ncon
    trials = find(trial_use{s} == Ntrial * (c-1) + 1):find(trial_use{s} == Ntrial * c);
    if ~isempty(trials)
        t = 1;
        for tr = trials
            EP_T = et.Epoch{tr};
            EP_Te.dat = EP_T.dat(:, :, 1:4 * Nrepeat);
            EP_Te.lat = EP_T.lat(:, 1:4 * Nrepeat);

            param_temp.trD.mode = 'testing';
            EP_test = Epoch_condition(EP_Te, param_temp);

            wind = param_temp.Baseline+1:param_temp.Totalepoc;
            EP_test.nar = preprocess_epochs(EP_test.nar, wind);
            test_temp = reshape_epochs(EP_test.nar);

            [rho1, rho2] = compute_rhos(test_temp, param_temp, Model_Init, Windowlist_new, Windowlength, Ntime_new);

            [score_1,score_2,pred_1, pred_2, predicted] = classify_rhos(Model_Init.model, rho1, rho2);

            if c > 1
                [c1, c2] = get_max_predictions(pred_1, pred_2);

                [upcheck1, ~, ~, ~, D1] = checkupdate(score_1, [], [], c1, []);
                [upcheck2, ~, ~, ~, D2] = checkupdate(score_2, [], [], c2, []);

                if should_update_model(adaptmode, fusionmode, upcheck1, upcheck2)
                    [param_temp, Model_Init.model] = update_model(EP_test, param_temp, Ntime_new, fusionmode, c1, c2, upcheck1, upcheck2, Model_Init.model);
                    fprintf('Block %d updated\n',tr);
                else
                    fprintf('Block %d NOT updated\n',tr);
                end


                [DSPs, Lambdas, LDA, UpDistance1, UpDistance2, Updated1, Updated2] = store_update_results(DSPs, Lambdas, LDA, UpDistance1, UpDistance2, Updated1, Updated2, t, c, s, param_temp, D1, D2, upcheck1, upcheck2, Ntime_new);
            end

            [scores1, scores2, outputs1, outputs2, outputs, answers] = store_results(scores1, scores2, outputs1, outputs2, outputs, answers, rho1, rho2, pred_1, pred_2, predicted, et, t, c, s, tr);

            t = t + 1;
        end
    end
end
end

function nar = preprocess_epochs(nar, wind)
nar = nar(wind, :, :, :);
nar = nar - mean(nar, 1, 'omitnan');
end

function reshaped = reshape_epochs(nar)
reshaped = permute(nar, [2, 1, 3, 4]);
reshaped = reshape(reshaped, size(reshaped, 1), size(reshaped, 2), []);
end

function [rho1, rho2] = compute_rhos(test_temp, param_temp, Model_Init, Windowlist_new, Windowlength, Ntime_new)
rho1 = NaN(size(test_temp, 3), Ntime_new);
rho2 = NaN(size(test_temp, 3), Ntime_new);

for ti = 1:Ntime_new
    windowind = Windowlist_new(ti):5:Windowlist_new(ti) + Windowlength - 1;

    W1 = param_temp.DSP(ti).W(:, 1 : param_temp.DSP(ti).nf);
    W2 = Model_Init.DSP(ti).W(:, 1 : param_temp.DSP(ti).nf);

    for tt = 1:size(test_temp, 3)
        Y_1 = W1' * test_temp(:, windowind, tt);
        Y_2 = W2' * test_temp(:, windowind, tt);

        [sc_11] = pred_mj(param_temp.trD.mdl(ti), Y_1(:)');
        [sc_112] = pred_mj(Model_Init.LDA(ti), Y_2(:)');

        rho1(tt, ti) = sc_11(:, end);
        rho2(tt, ti) = sc_112(:, end);
    end
end
end

function [score_1,score_2,pred_1, pred_2, predicted] = classify_rhos(model, rho1, rho2)
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

function [param, model] = update_model(EP_update, param, Ntime_new, fusionmode, c1, c2, upcheck1, upcheck2, model)
if strcmp(fusionmode, 'fusion')
    if upcheck1
        uplabel = c1;
    elseif upcheck2
        uplabel = c2;
    end
else
    uplabel = c1;
end


Feat = [];
yynew = [];
for ti = 1:Ntime_new
    param_temp = param;
    param_temp.DSP  = param.DSP(ti);
    param_temp = updateDSP(EP_update, uplabel, param_temp);
    param.DSP(ti) = param_temp.DSP;

    Feat = FeatureExt_DSP(EP_update, param_temp);
    LDA = updateSingleClassfier(Feat, param.trD.mdl(ti), []);
    param.trD.mdl(ti) = LDA;
    [sc_2] = pred_mj(param.trD.mdl(ti), Feat);
    yynew(:, ti) = sc_2(:, end);
end

%     if size(Xnew, 1) < Ntime * 50
%         label_temp = zeros(size(Feat{1}, 1), 1);
%         label_temp = reshape(label_temp, Nrepeat, 4);
%         label_temp(:, uplabel) = 1;
%         label = label_temp(:);
%         Xnew = [Xnew; yynew];
%         Ynew = [Ynew; label];
%         model = fitglm(Xnew, Ynew, 'linear', 'Distribution', 'binomial');
%     end
end

function [DSPs, Lambdas, LDA, UpDistance1, UpDistance2, Updated1, Updated2] = store_update_results(DSPs, Lambdas, LDA, UpDistance1, UpDistance2, Updated1, Updated2, t, c, s, param_temp, D1, D2, upcheck1, upcheck2, Ntime_new)
dsptemp = [];
for ti = 1:Ntime_new
    dsptemp(:, :, ti) = param_temp.DSP(ti).W;
    lambdas(ti) = param_temp.DSP(ti).lambda;
end
DSPs{t, c, s} = dsptemp;
Lambdas(:,t, c, s) = lambdas;
LDA{t, c, s} = param_temp.trD.mdl;
UpDistance1(:, t, c, s) = D1;
UpDistance2(:, t, c, s) = D2;
Updated1(t, c, s) = upcheck1;
Updated2(t, c, s) = upcheck2;
end

function [scores1, scores2, outputs1, outputs2, outputs, answers] = store_results(scores1, scores2, outputs1, outputs2, outputs, answers, rho1, rho2, pred_1, pred_2, predicted, et, t, c, s, tr)
scores1(:, :, t, c, s) = rho1;
scores2(:, :, t, c, s) = rho2;
outputs1(:, t, c, s) = pred_1;
outputs2(:, t, c, s) = pred_2;
outputs(t, c, s) = predicted;
answers(t, c, s) = et.Epoch{tr}.target;
end

function save_results(fpathanalysis, adaptmode, fusionmode, scores1, scores2, outputs1, outputs2, outputs, answers, DSPs, BadChs, LDA, Updated1, Updated2, UpDistance1, UpDistance2, Lambdas, mdl_init)
cd(fpathanalysis);
mkdir('NewFunc_adaptiveMWHDPA_LDA_MultiOnly_NoThre_VarNf_VarUC_0.06');
save(['NewFunc_adaptiveMWHDPA_LDA_MultiOnly_NoThre_VarNf_VarUC_0.06\', adaptmode, '_', fusionmode, '.mat'], ...
    'scores1', 'scores2', 'outputs1', 'outputs2', 'outputs', 'answers', ...
    'DSPs', 'BadChs', 'LDA', 'Updated1', 'Updated2', 'UpDistance1', 'UpDistance2', 'Lambdas', 'mdl_init');
end
