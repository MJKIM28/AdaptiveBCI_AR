function [hitLabel,param]= Classification(Feature,label,param)
global gamma
if strcmp(param.trD.mode,'training')
    try
        gamma = 1/(size(Feature,2)*var(Feature(:)));
        
        param.trD.mdl = fitcdiscr(Feature, label);%,'KernelFunction','mysigmoid');
%         param.trD.mdl = fitPosterior(param.trD.mdl);
        
        param.trD.gamma = gamma;
        
        param.trD.mode      = 'testing';
        
        
%         [~,score] = predict(param.trD.mdl,Feature);

%         
%         param.trD.SV = param.trD.mdl.X(param.trD.mdl.IsSupportVector,:);
%         param.trD.Y_SV = label(param.trD.mdl.IsSupportVector,:);
%         param.trD.score_SV = score(param.trD.mdl.IsSupportVector,:);
%         param.trD.Alpha = param.trD.mdl.Alpha;
%         param.trD.Bias = param.trD.mdl.Bias;
%         param.trD.score = score;

        hitLabel = [];
        
        Cvind =      crossvalind('Kfold',label,5);
        fprintf('\n');
        for r = 1:5
            fprintf('%d',r)
            trind = Cvind ~=r ;
            feat = Feature(trind,:);
            gamma = 1/(size(feat,2)*var(feat(:)));
            mdl=  fitcdiscr(feat,label(trind));%,'KernelFunction','mysigmoid');
            teind = Cvind == r;
            [result, score] = predict(mdl,Feature(teind,:));
            acc(r) = length(find(result == label(teind)))/length(result);
        end
        fprintf('\n');
        CVAcc = mean(acc(r));
        fprintf('CV accuracy: %.2f\n',CVAcc);
        param.CVAcc = CVAcc;
        
    catch
        fprintf('Training failed..!\n')
        hitLabel = -1;
    end
else
    try
        
        %         Feature = (Feature - repmat( param.trD.featparam(1), size(Feature,1),1))./repmat(param.trD.featparam(2), size(Feature,1),1);
        gamma = param.trD.gamma;
%         gamma = 1/(size(Feature,2)*var(Feature(:)));
        [C,sc] = predict(param.trD.mdl, Feature);
        if isfield(param,'DSP') || isfield(param.DTP,'W')
            if isfield(param.DSP,'W') || isfield(param.DTP,'W')
                
%                 for tr = 1:size(Feature,1)
%                     % check if ERP amplitude was 0 (ignored trial, only DSP w0 values) 
%                    if sum(Feature(tr,:) == reshape(param.DSP.w0(1:param.DSP.nf,:),1,size(Feature,2))) == size(Feature,2)
%                        sc(tr,:) = 0; 
%                    end
%                 end

                Nt = size(Feature,1)/(param.repeat*param.NumStims);
                tescore_bl = reshape( sc(:,end),param.repeat,param.NumStims,Nt);
                tescore_bl_sum = sum(tescore_bl,1);
                [~,hit_id] = max(tescore_bl_sum,[],2);
                hit_id = squeeze(hit_id);
                param.trD.score = tescore_bl;
            end
        else
            [~,hit_id] = max(sc(:,end));
        end
        hitLabel = param.Stims(hit_id);
        
        %% for DS
        % If iteration end but block not end
        % accumulate classifier output
%         if (~param.switch_on) || isfield(param,'switch_on_iter')
%             param.decoder.data = [param.decoder.data; sc(:,end)];
%         end
        
    catch
        hitLabel=-1;
    end
end