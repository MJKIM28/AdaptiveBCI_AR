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

    EP = EpochData(sig,trig,param);
    [Feature, label, param] = FeatureExtraction(EP, param);

    switch param.trD.mode
        case 'training'
            %-- train classifier
            [C,param]               = Classification(Feature,label, param);
        case 'testing'
            %-- initial output
            [C,param] = Classification(Feature,label,param); % use updated Classifier
            score = param.trD.score;
            output = getoutput(score); %SVM kernel

            %-- check update
            [up,Posterior,D] = checkupdate(score, [] ,param.trD.threshold,output,[]);
            param.update{param.Numtrial}.updated = up;
            param.update{param.Numtrial}.Posterior = Posterior;
            param.update{param.Numtrial}.Dist = D;

            %-- update
            %-- DSP
            if ~isempty(up)%%&& con ~= 1
                EP_1block = Epoch_condition(EP,param);

                %-- update DSP
                if length(up) > 1
                    EP_update = EP_1block;
                    EP_update.nar = EP_update.nar(:,:,up,:);

                else
                    EP_update = EP_1block;
                end

                
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

                param.update{param.Numtrial}.DSP = param.DSP;
                param.update{param.Numtrial}.mdl = param.trD;
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