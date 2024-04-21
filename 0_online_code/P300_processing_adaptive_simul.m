function [C, param]   = P300_processing_adaptive(sig, trig, param)
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


     EP                      = Epoching_singletr_simul(sig, trig, param);
    param.Target = EP.target;
            
    [Feature, label, param] = FeatureExtraction(EP, param);

    switch param.trD.mode
        case 'training'
            %-- train classifier
            [C,param]               = Classification(Feature,label, param);
            param.trD.mdl_init = param.trD.mdl;
            param.trD.mdl_adapt = param.trD.mdl;
            param.trD.feature = Feature;
            param.trD.label = label;
            param.update(1).DSP_init = param.DSP;

            param.DSP.init = param.DSP;
        case 'testing'

            if strcmp(param.trD.ADmode,'fixed')
                param.trD.mdl = param.trD.mdl_init;
            [C,param] = Classification(Feature,label,param); % use updated Classifier
                     fprintf('>> Use fixed\n') 
      

            elseif strcmp(param.trD.ADmode,'adaptive')
                param.trD.mdl = param.trD.mdl_adapt;
                 %-- initial output
            [C,param] = Classification(Feature,label,param); % use updated Classifier
            score = param.trD.score;
            output = getoutput(score); %SVM kernel

            
                %-- check update
                [upcheck,up,Posterior,trs,D] = checkupdate(score, [] ,param.trD.threshold,output,[]);
               
                param.update(param.Numtrial).classprob = score;
                param.update(param.Numtrial).upcheck = upcheck;
                param.update(param.Numtrial).upids = up;
                param.update(param.Numtrial).trial_candidate = trs;
                param.update(param.Numtrial).Posterior = Posterior;
                param.update(param.Numtrial).Dist = D;
                param.update(param.Numtrial).prediction = C;

                %-- update
                %-- DSP
                if upcheck && ~isempty(up)
                    EP_1block = Epoch_condition(EP,param);

                    %-- update DSP
                    EP_update = EP_1block;
                    EP_update.nar = EP_update.nar(:,:,up,:);

                    param  = updateDSP(EP_update,output,param);

                    %-- apply updated DSP to feat

                    Repeat = param.repeat;
                    param.repeat = size(EP_update.nar,3);
                    Feat = FeatureExt_DSP(EP_update,param);
                    [C_new,param] = Classification(Feat,[],param);

                    label_temp = -ones(size(Feat,1),1);
                    label_temp = reshape(label_temp,param.repeat,param.NumStims);
                    label_temp(:,C_new) = 1;
                    label = label_temp(:);

                    %-- re-calibration
                    param.trD.feature = [param.trD.feature; Feat];
                    param.trD.label = [param.trD.label; label];

                    param.trD.mode = 'training';
                    [~,param] = Classification(param.trD.feature,param.trD.label,param);

                    param.repeat = Repeat;

                    fprintf('>> Updated\n')


                    param.update(param.Numtrial).DSP = param.DSP;
                    param.update(param.Numtrial).mdl = param.trD;
                    param.update(param.Numtrial).prediction_new = C_new;
                    param.trD.mdl_adapt = param.trD.mdl;
                else
                    fprintf('>> Not updated\n')
                    param.update(param.Numtrial).prediction_new = NaN;
                end
            
            end

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