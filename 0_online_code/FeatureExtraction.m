function [Feature, label, param] = FeatureExtraction(EP, param)
global gamma
try
    switch param.FeatureType
        case 'waveform' % type 1
            if strcmp(param.EpocType,'singletr') || strcmp(param.EpocType,'singletr_tsb')
                [EP_b, param] = FeatExt_basic_prep(EP,param);
            elseif strcmp(param.EpocType,'meantr') || strcmp(param.EpocType,'meantr_tsb')
                EP_b = EP;
            end
            [Feature,label,param] = FeatureExt_basic(EP_b,param);

        case 'DSP' % type 2
            if strcmp(param.EpocType,'singletr') || strcmp(param.EpocType,'singletr_tsb')
                EP = Epoch_condition(EP,param);

                if strcmp(param.trD.mode,'training')
                    tic
                    nf_param = 1:ceil(param.NumCh/2);
                    %--  DSP param select
                    CV = [];

                    n = 1;
                    for nf = nf_param
                        param.DSP.nf = nf;
                        if isfield(param.DSP,'W')
                            param.DSP = rmfield(param.DSP,'W');
                        end
                        [Feature, label, param] = FeatureExt_DSP(EP,param);
                        Cvind =      crossvalind('Kfold',label,5);
                        acc=[];
                        for r = 1:5
                            trind = Cvind ~=r ;
                            feat = Feature(trind,:);
                            gamma = 1/(size(feat,2)*var(feat(:)));
                            mdl=  fitcsvm(feat,label(trind),'KernelFunction','mysigmoid');
                            teind = Cvind == r;
                            [result, score] = predict(mdl,Feature(teind,:));
                            acc(r) = length(find(result == label(teind)))/length(result);
                        end
                        fprintf('.');
                        CVAcc = mean(acc(r));
                        CV(n) = CVAcc;
                        n = n+1;
                    end
                    fprintf('\n');
                    toc
                    [cv_max1,nf_sel_id] = max(CV);
                    param.DSP.nf = nf_param(nf_sel_id);
                    if isfield(param.DSP,'W')
                        param.DSP = rmfield(param.DSP,'W');
                    end
                end
                [Feature, label, param] = FeatureExt_DSP(EP,param);

            end
    end

    % type 3
    % [Feature, label, param] = FeatureExt_xDAWN(EP,param);

catch
    keyboard
end