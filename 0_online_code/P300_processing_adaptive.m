function [C, param]   = P300_processing_adaptive(sig, trig, param)
global DS
% updated @ 20240311

try
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

    switch param.trD.mdl
        case 'training'
            %-- train classifier
            [C,param]               = Classification(Feature,label, param);
        case 'testing'
            %-- initial output
            [C,param] = Classification_new(Feature,label,param); % use updated Classifier
            score = param.trD.score;
            output = getoutput(score); %SVM kernel

            %-- check update
            [up,Posterior] = checkupdate(score, [] ,threshold,output);


            %-- update
            %-- DSP
            if up
                EP_1block = Epoch_condition(EP,paramL);
                param  = updateDSP(EP_1block,output,param);
                %-- apply updated DSP to feat
                Feat = FeatureExt_DSP(EP_1block,param);
                %-- SVM
                param = updateClassifier(Feat,output,param);
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
catch
    keyboard
end