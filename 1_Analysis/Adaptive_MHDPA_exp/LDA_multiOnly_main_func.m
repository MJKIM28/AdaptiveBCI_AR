function LDA_main(adaptmode, fusionmode, k)
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
    Windowlist = fix(1:0.05 * Fs:0.6 * Fs - Windowlength);
    Ntime = length(Windowlist);
    chanloc = readlocs('D:\[1]EEGBCI\[2]Research\Code&Algorithm\BP_channelLocs.locs');
    trial_use = load_trials();

    % Initialize storage variables
    [scores1, scores2, outputs1, outputs2, outputs, answers, DSPs, BadChs, Lambdas, LastDecision, UpDistance1, UpDistance2, Updated1, Updated2, LDA] = initialize_storage(Ntime, Nrepeat, Ncon, Nsub);

    for s = 1:Nsub
        SubName = SubNames{s};
        fprintf('%s\n', SubName);

        [e, et, param, chlist] = load_subject_data(fpath, SubName, chanloc);
        [param, EP_train_tr, EP_train_val, val_temp, Label_val, Nval] = prepare_training_data(e, param, Nrepeat, Fs);

        % Training Part with k-fold cross-validation
        [Model_Init, param_temp, Ntime_new, Windowlist_new] = train_model_kfold(param, EP_train_tr, val_temp, Label_val, chlist, Windowlist, Windowlength, Ntime, k);

        % Testing Part
        [scores1, scores2, outputs1, outputs2, outputs, answers, DSPs, Lambdas, LDA, UpDistance1, UpDistance2, Updated1, Updated2] = ...
            test_blocks(adaptmode, fusionmode, Ncon, trial_use, s, Nrepeat, et, param_temp, Model_Init, Ntime_new, Windowlist_new, ...
                        scores1, scores2, outputs1, outputs2, outputs, answers, DSPs, Lambdas, LDA, UpDistance1, UpDistance2, Updated1, Updated2);
    end

    save_results(fpathanalysis, adaptmode, fusionmode, scores1, scores2, outputs1, outputs2, outputs, answers, DSPs, BadChs, LDA, Updated1, Updated2, UpDistance1, UpDistance2, Lambdas, LastDecision);
end

function trial_use = load_trials()
    load('E:\1_연구\분석코드\2024_adaptive_BCI\AdaptiveBCI_AR\1_Analysis\TrialUsed2.mat');
    load('E:\1_연구\분석코드\2024_adaptive_BCI\AdaptiveBCI_AR\1_Analysis\TrialUsed.mat');
    trial_use1 = Trials;
    trial_use1 = trial_use1(9:end);
    trial_use2 = Trials2;
    trial_use = [trial_use1; trial_use2];
end

function [scores1, scores2, outputs1, outputs2, outputs, answers, DSPs, BadChs, Lambdas, LastDecision, UpDistance1, UpDistance2, Updated1, Updated2, LDA] = initialize_storage(Ntime, Nrepeat, Ncon, Nsub)
    scores1 = NaN(40, Ntime, Nrepeat, Ncon, Nsub);
    scores2 = NaN(40, Ntime, Nrepeat, Ncon, Nsub);
    outputs1 = NaN(4, Nrepeat, Ncon, Nsub);
    outputs2 = NaN(4, Nrepeat, Ncon, Nsub);
    outputs = NaN(Nrepeat, Ncon, Nsub);
    answers = NaN(Nrepeat, Ncon, Nsub);
    DSPs = [];
    BadChs = cell(Nsub, 1);
    Lambdas = cell(Nrepeat, Ncon, Nsub);
    LDA = cell(Nrepeat, Ncon, Nsub);
    Updated1 = NaN(Nrepeat, Ncon, Nsub);
    Updated2 = NaN(Nrepeat, Ncon, Nsub);
    UpDistance1 = NaN(Nrepeat, Ncon, Nsub);
    UpDistance2 = NaN(Nrepeat, Ncon, Nsub);
    LastDecision = cell(Nsub, 1);
end

function [e, et, param, chlist] = load_subject_data(fpath, SubName, chanloc)
    e = load([fpath, '\Train\Epoch\', SubName]);
    et = load([fpath, '\Test\Epoch\', SubName]);
    p = load([fpath, '\Dat_', SubName, '\param.mat']);
    param = p.param;
    param.DSP = rmfield(param.DSP, 'W');
    param.chanlocIn = chanloc;
    param.chanlocIn(param.badch) = [];
    chlist = setdiff(1:length(chanloc), param.badch);
end

function [param, EP_train_tr, EP_train_val, val_temp, Label_val, Nval] = prepare_training_data(e, param, Nrepeat, Fs)
    param.trD.mode = 'training';
    param.repeat = Nrepeat;
    EP = e.Epoch;
    [EP_train_tr, EP_train_val] = segment_epochs(EP, Nrepeat, Fs);
    [val_temp, Label_val, Nval] = prepare_validation_data(EP_train_val, Fs);
end

function [EP_train_tr, EP_train_val] = segment_epochs(EP, Nrepeat, Fs)
    trid = 1:10;
    valid = 11:15;
   

end

function [val_temp, Label_val, Nval] = prepare_validation_data(EP_train_val, Fs)
    val_temp = cat(3, EP_train_val.tar, EP_train_val.nar);
    val_temp = permute(val_temp, [2, 1, 3]);
    Nval = size(val_temp, 3);
    Label_val = [ones(size(EP_train_val.tar, 3), 1); zeros(size(EP_train_val.nar, 3), 1)];
end

function [Model_Init, param_temp, Ntime_new, Windowlist_new] = train_model_kfold(param, EP_train_tr, val_temp, Label_val, chlist, Windowlist, Windowlength, Ntime, k)
    folds = crossvalind('Kfold', size(EP_train_tr.target, 2), k);
    Acc_val = zeros(fix(length(chlist) / 2), Ntime, k);

    for fold = 1:k
        train_idx = folds ~= fold;
        val_idx = folds == fold;

        EP_train_tr_fold = select_epochs(EP_train_tr, train_idx,param.repeat);
        EP_train_val_fold = select_epochs(EP_train_tr, val_idx,param.repeat);

        val_temp_fold = prepare_validation_data(EP_train_val_fold, Fs);

        [Acc_val(:, :, fold), Nfs_fold] = validate_NF(param, EP_train_tr_fold, val_temp_fold, Label_val(val_idx), chlist, Windowlist, Windowlength, Ntime);
    end

    Acc_val = mean(Acc_val, 3);
    [~, Nfs] = max(Acc_val);

    [Model_Init, param_temp, Ntime_new, Windowlist_new] = initialize_model(param, EP_train_tr, val_temp, Label_val, Nfs, Windowlist, Windowlength, Ntime);
end

function EP_selected = select_epochs(EP, idx,Nrepeat)
 EP_selected.target = EP.target(idx);
    tridall = bsxfun(@plus, (idx-1) * Nrepeat, [1:Nrepeat]');
    tridall2 = bsxfun(@plus, (idx-1) * Nrepeat * 3, [1:Nrepeat * 3]');
    EP_selected.tar = EP.tar(0.2 * Fs + 1:0.8 * Fs, :, tridall(:));
    EP_selected.nar = EP.nar(0.2 * Fs + 1:0.8 * Fs, :, tridall2(:));


    EP_selected.tar = EP_selected.tar - mean(EP_selected.tar, 1, 'omitnan');
    EP_selected.nar = EP_selected.nar - mean(EP_selected.nar, 1, 'omitnan');

  
end

function [Acc_val, Nfs] = validate_NF(param, EP_train_tr, val_temp, Label_val, chlist, Windowlist, Windowlength, Ntime)
    Acc_val = [];
    for Nf = 1:fix(length(chlist) / 2)
        param.DSP.nf = Nf;
        yy = [];
        param_temp = [];
        for ti = 1:Ntime
            windowind = Windowlist(ti):5:Windowlist(ti) + Windowlength - 1;
            param.DSP.window = windowind;
            [~, ~, param_temp{ti}] = FeatureExt_DSP(EP_train_tr, param);

            [Template_Tar, Template_Nar] = filter_templates(param_temp{ti}, EP_train_tr, windowind, Nf);
            X = [reshape(Template_Tar, [], size(Template_Tar, 3))'; reshape(Template_Nar, [], size(Template_Nar, 3))'];
            Y = [ones(size(Template_Tar, 3), 1); zeros(size(Template_Nar, 3), 1)];

            mdl_mj{ti} = getLDA(X, Y);
            X_ = extract_validation_features(param_temp{ti}, val_temp, windowind);
            [sc_2] = pred_mj(mdl_mj{ti}, reshape(X_, [], size(X_, 3))');
            yy(:, ti) = sc_2(:, end);
        end
        Acc_val(Nf, :) = sum(yy > 0.5 == Label_val) / length(Label_val);
    end
    [~, Nfs] = max(Acc_val);
end

function [Template_Tar, Template_Nar] = filter_templates(param_temp, EP_train_tr, windowind, Nf)
    Template_Tar = [];
    Template_Nar = [];
    for tr = 1:size(EP_train_tr.tar, 3)
        Template_Tar(:, :, tr) = param_temp.DSP.W(:, 1:Nf)' * EP_train_tr.tar(windowind, :, tr)';
    end
    for tr = 1:size(EP_train_tr.nar, 3)
        Template_Nar(:, :, tr) = param_temp.DSP.W(:, 1:Nf)' * EP_train_tr.nar(windowind, :, tr)';
    end
end

function X_ = extract_validation_features(param_temp, val_temp, windowind)
    X_ = [];
    for tr = 1:size(val_temp, 3)
        X_(:, :, tr) = param_temp.DSP.W(:, 1:param_temp.DSP.nf)' * val_temp(:, windowind, tr);
    end
end

function [Model_Init, param_temp, Ntime_new, Windowlist_new] = initialize_model(param, EP_train_tr, val_temp, Label_val, Nfs, Windowlist, Windowlength, Ntime)
    yy = [];
    param_temp = [];
    for ti = 1:Ntime
        windowind = Windowlist(ti):5:Windowlist(ti) + Windowlength - 1;
        param.DSP.window = windowind;
        param.DSP.nf = Nfs(ti);

        [~, ~, param_temp{ti}] = FeatureExt_DSP(EP_train_tr, param);
        [Template_Tar, Template_Nar] = filter_templates(param_temp{ti}, EP_train_tr, windowind, Nfs(ti));
        X = [reshape(Template_Tar, [], size(Template_Tar, 3))'; reshape(Template_Nar, [], size(Template_Nar, 3))'];
        Y = [ones(size(Template_Tar, 3), 1); zeros(size(Template_Nar, 3), 1)];

        mdl_mj{ti} = getLDA(X, Y);
        X_ = extract_validation_features(param_temp{ti}, val_temp, windowind);
        [sc_2] = pred_mj(mdl_mj{ti}, reshape(X_, [], size(X_, 3))');
        yy(:, ti) = sc_2(:, end);

        param_temp{ti}.DSP.lambda = NaN;
        param_temp{ti}.trD.mode = 'testing';
        mdl_mj{ti}.lambda = NaN;
        param_temp{ti}.DSP.CosSim = NaN(length(chlist), 1);
        mdl_mj{ti}.CosSim = NaN;
    end

    X_init = yy;
    Y_init = Label_val;
    model = fitglm(X_init, Y_init, 'linear', 'Distribution', 'binomial');

    Model_Init.DSP = param_temp;
    Model_Init.LDA = mdl_mj;
    Model_Init.model = model;
    Windowlist_new = Windowlist;
    Ntime_new = Ntime;
end

function [scores1, scores2, outputs1, outputs2, outputs, answers, DSPs, Lambdas, LDA, UpDistance1, UpDistance2, Updated1, Updated2] = test_blocks(adaptmode, fusionmode, Ncon, trial_use, s, Nrepeat, et, param_temp, Model_Init, Ntime_new, Windowlist_new, scores1, scores2, outputs1, outputs2, outputs, answers, DSPs, Lambdas, LDA, UpDistance1, UpDistance2, Updated1, Updated2)
    for c = 1:Ncon
        trials = find(trial_use{s} == Nrepeat * (c-1) + 1):find(trial_use{s} == Nrepeat * c);
        if ~isempty(trials)
            t = 1;
            for tr = trials
                EP_T = et.Epoch{tr};
                EP_Te.dat = EP_T.dat(:, :, 1:4 * Nrepeat);
                EP_Te.lat = EP_T.lat(:, 1:4 * Nrepeat);

                param.trD.mode = 'testing';
                EP_test = Epoch_condition(EP_Te, param);

                EP_test.nar = preprocess_epochs(EP_test.nar, Fs);
                test_temp = reshape_epochs(EP_test.nar);

                [rho1, rho2] = compute_rhos(test_temp, param_temp, Model_Init, Windowlist_new, Windowlength, Ntime_new);

                [pred_1, pred_2, predicted] = classify_rhos(Model_Init.model, rho1, rho2);

                if c > 1
                    [c1, c2] = get_max_predictions(pred_1, pred_2);

                    [upcheck1, D1] = check_update(pred_1, c1);
                    [upcheck2, D2] = check_update(pred_2, c2);

                    if should_update_model(adaptmode, fusionmode, upcheck1, upcheck2)
                        [param_temp, mdl_mj, Model_Init.model] = update_model(EP_test, param_temp, Ntime_new, fusionmode, c1, c2, upcheck1, upcheck2, Nrepeat, Windowlength, Windowlist_new, mdl_mj, Model_Init.model);
                    end

                    [DSPs, Lambdas, LDA, UpDistance1, UpDistance2, Updated1, Updated2] = store_update_results(DSPs, Lambdas, LDA, UpDistance1, UpDistance2, Updated1, Updated2, t, c, s, param_temp, mdl_mj, D1, D2, upcheck1, upcheck2, Ntime_new);
                end

                [scores1, scores2, outputs1, outputs2, outputs, answers] = store_results(scores1, scores2, outputs1, outputs2, outputs, answers, rho1, rho2, pred_1, pred_2, predicted, et, t, c, s, tr);

                t = t + 1;
            end
        end
    end
end

function nar = preprocess_epochs(nar, Fs)
    nar = nar(0.2 * Fs + 1 : 0.8 * Fs, :, :, :);
    nar = nar - mean(nar, 1, 'omitnan');
end

function reshaped = reshape_epochs(nar)
    reshaped = permute(nar, [2, 1, 3, 4]);
    reshaped = reshape(reshaped, size(reshaped, 1), size(reshaped, 2), []);
end

function [rho1, rho2] = compute_rhos(test_temp, param_temp, Model_Init, Windowlist_new, Windowlength, Ntime_new)
    rho1 = [];
    rho2 = [];
    for tt = 1:size(test_temp, 3)
        for ti = 1:Ntime_new
            windowind = Windowlist_new(ti):5:Windowlist_new(ti) + Windowlength - 1;

            Y_1 = param_temp{ti}.DSP.W(:, 1 : param_temp{ti}.DSP.nf)' * test_temp(:, windowind, tt);
            [sc_11] = pred_mj(mdl_mj{ti}, Y_1(:)');
            rho1(tt, ti) = sc_11(:, end);

            Y_2 = Model_Init.DSP{ti}.DSP.W(:, 1 : param_temp{ti}.DSP.nf)' * test_temp(:, windowind, tt);
            [sc_112] = pred_mj(Model_Init.LDA{ti}, Y_2(:)');
            rho2(tt, ti) = sc_112(:, end);
        end
    end
end

function [pred_1, pred_2, predicted] = classify_rhos(model, rho1, rho2)
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

function [upcheck, D] = check_update(pred, c)
    [upcheck, ~, ~, ~, D] = checkupdate(pred, [], [], c, []);
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

function [param_temp, mdl_mj, model] = update_model(EP_update, param_temp, Ntime_new, fusionmode, c1, c2, upcheck1, upcheck2, Nrepeat, Windowlength, Windowlist_new, mdl_mj, model)
    if strcmp(fusionmode, 'fusion')
        if upcheck1
            uplabel = c1;
        elseif upcheck2
            uplabel = c2;
        end
    else
        uplabel = c1;
    end
    fprintf('Block %d updated\n', tr);

    Feat = [];
    yynew = [];
    for ti = 1:Ntime_new
        param_temp{ti} = updateDSP(EP_update, uplabel, param_temp{ti});
        Feat{ti} = FeatureExt_DSP(EP_update, param_temp{ti});
        mdl_mj{ti} = updateSingleClassfier(Feat{ti}, mdl_mj{ti}, []);
        [sc_2] = pred_mj(mdl_mj{ti}, Feat{ti});
        yynew(:, ti) = sc_2(:, end);
    end

    if size(Xnew, 1) < Ntime * 50
        label_temp = zeros(size(Feat{1}, 1), 1);
        label_temp = reshape(label_temp, Nrepeat, 4);
        label_temp(:, uplabel) = 1;
        label = label_temp(:);
        Xnew = [Xnew; yynew];
        Ynew = [Ynew; label];
        model = fitglm(Xnew, Ynew, 'linear', 'Distribution', 'binomial');
    end
end

function [DSPs, Lambdas, LDA, UpDistance1, UpDistance2, Updated1, Updated2] = store_update_results(DSPs, Lambdas, LDA, UpDistance1, UpDistance2, Updated1, Updated2, t, c, s, param_temp, mdl_mj, D1, D2, upcheck1, upcheck2, Ntime_new)
    dsptemp = [];
    for ti = 1:Ntime_new
        dsptemp(:, :, ti) = param_temp{ti}.DSP.W;
        lambdas(ti) = param_temp{ti}.DSP.lambda;
    end
    DSPs{t, c, s} = dsptemp;
    Lambdas{t, c, s} = lambdas;
    LDA{t, c, s} = mdl_mj;
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

function save_results(fpathanalysis, adaptmode, fusionmode, scores1, scores2, outputs1, outputs2, outputs, answers, DSPs, BadChs, LDA, Updated1, Updated2, UpDistance1, UpDistance2, Lambdas, LastDecision)
    cd(fpathanalysis);
    mkdir('adaptiveMWHDPA_LDA_MultiOnly_NoThre_VarNf_VarUC_updateLRmax50');
    save(['adaptiveMWHDPA_LDA_MultiOnly_NoThre_VarNf_VarUC_updateLRmax50\', adaptmode, '_', fusionmode, '.mat'], ...
        'scores1', 'scores2', 'outputs1', 'outputs2', 'outputs', 'answers', ...
        'DSPs', 'BadChs', 'LDA', 'Updated1', 'Updated2', 'UpDistance1', 'UpDistance2', 'Lambdas', 'LastDecision');
end
