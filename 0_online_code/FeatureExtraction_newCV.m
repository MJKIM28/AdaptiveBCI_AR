function [Feature, label, param] = FeatureExtraction(EP, param)
global gamma
Nfold = 5;
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
                nrepeat = param.repeat;
                if strcmp(param.trD.mode,'training')

                    tic
                    nf_param = 1:ceil(param.NumCh/2);
                    %% %--  DSP param select
                    CV = NaN(1,length(nf_param));

                    
                    fprintf('CV\n')
                    parfor nf = nf_param
                       Cvind =      crossvalind('Kfold',length(EP.target),Nfold);
                        acc=[];
                        
                        for r = 1:Nfold
                            
                            param_temp = param;
                            param_temp.DSP.nf = nf;
                            trind = find(Cvind ~=r) ;
                            EP_train = [];
                            EP_train.target = EP.target(trind);
                            trialind_tar = bsxfun(@plus,(trind-1)*nrepeat, 1:nrepeat)';
                            trialind_nar = bsxfun(@plus,(trind-1)*nrepeat*(param_temp.NumStims-1), 1:nrepeat*(param_temp.NumStims-1))';
                            EP_train.tar = EP.tar(:,:,trialind_tar(:));
                            EP_train.nar = EP.nar(:,:,trialind_nar(:));
                            param_temp.trD.mode = 'training';
                            if isfield(param_temp.DSP,'W')
                                param_temp.DSP = rmfield(param_temp.DSP,'W');
                            end
                            [feat, label, param_temp] = FeatureExt_DSP(EP_train,param_temp);


                            [mdl,~,~] = classification_inCV(feat,label,[],'train');

                            teind = find(Cvind == r);
                            EP_test = [];
                            EP_test.target = EP.target(teind);
                            trialind_te_tar = bsxfun(@plus,(teind-1)*nrepeat, 1:nrepeat)';
                            trialind_te_nar = bsxfun(@plus,(teind-1)*nrepeat*(param_temp.NumStims-1), 1:nrepeat*(param_temp.NumStims-1))';
                            EP_test.tar = EP.tar(:,:,trialind_te_tar(:));
                            EP_test.nar = EP.nar(:,:,trialind_te_nar(:));
                            
                            [feat_te, label_te, param_temp] = FeatureExt_DSP(EP_test,param_temp);
                              

                            [~,C,sc] = classification_inCV(feat_te,[],mdl,'test');
  
                            
                            acc(r) = length(find(C == label_te))/length(label_te);
                        end
                        fprintf('.');
                        CVAcc = mean(acc);
                        CV(nf) = CVAcc;
                    end
                    param.trD.mode = 'training';
                    fprintf('\n');
                    toc
                    [cv_max1,nf_sel_id] = max(CV);
                    fprintf('max CV acc: %.2f\n',cv_max1);
                    param.CVAcc_featext = CV;
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
end


function [mdl,C,sc] = classification_inCV(Feature,label,mdl,type)

global gamma

if strcmp(type,'train')
gamma = 1/(size(Feature,2)*var(Feature(:)));
 mdl=  fitcsvm(Feature,label);%,'KernelFunction','mysigmoid');
 C = 0; sc = 0;
else
    [C,sc] = predict(mdl,Feature);
end

% if strcmp(param.trD.mode,'training')
%     try
%         gamma = 1/(size(Feature,2)*var(Feature(:)));
%         
%         param.trD.mdl = fitcsvm(Feature, label,'KernelFunction','mysigmoid');
%         param.trD.gamma = gamma;
%         param.trD.mode      = 'testing';
%         hitLabel = [];
%         
%         
%     catch
%         fprintf('Training failed..!')
%         hitLabel = -1;
%     end
% else
%     try
%         
%                 Feature = (Feature - repmat( param.trD.featparam(1), size(Feature,1),1))./repmat(param.trD.featparam(2), size(Feature,1),1);
%         gamma = param.trD.gamma;
%         gamma = 1/(size(Feature,2)*var(Feature(:)));
%         [C,sc] = predict(param.trD.mdl, Feature);
%         if isfield(param,'DSP')
%             if isfield(param.DSP,'W')
%                 
%                 for tr = 1:size(Feature,1)
%                     check if ERP amplitude was 0 (ignored trial, only DSP w0 values) 
%                    if sum(Feature(tr,:) == reshape(param.DSP.w0(1:param.DSP.nf,:),1,size(Feature,2))) == size(Feature,2)
%                        sc(tr,:) = 0; 
%                    end
%                 end
%                 Nt = size(Feature,1)/(param.repeat*param.NumStims);
%                 tescore_bl = sum(reshape(sc(:,end),param.repeat,param.NumStims,Nt),1);
%                 [~,hit_id] = max(tescore_bl,[],2);
%                 hit_id = squeeze(hit_id);
%               
%             end
%         else
%             [~,hit_id] = max(sc(:,end));
%         end
%         hitLabel = param.Stims(hit_id);
        
        %% for DS
        % If iteration end but block not end
        % accumulate classifier output
%         if (~param.switch_on) || isfield(param,'switch_on_iter')
%             param.decoder.data = [param.decoder.data; sc(:,end)];
%         end
        
%     catch
%         hitLabel=-1;
%     end
end