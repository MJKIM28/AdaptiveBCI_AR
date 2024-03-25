function [hitLabel,param]= Classification_LR(Feature,label,param)
if strcmp(param.trD.mode,'training')
    try
  
%         param.trD.mdl = fitclinear(Feature,label,"Learner","logistic",'Regularization','lasso');
       
        if min(label) < 0
            label(label==-1)=0;
        end
        param.trD.mdl = fitglm(Feature,label,"Distribution","binomial",'Link','logit');
        [B,FitInfo] =lassoglm(Feature,label,"binomial", 'CV', 10);
        idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
        B0 = FitInfo.Intercept(idxLambdaMinDeviance); % Intercept
        coefficients = B(:, idxLambdaMinDeviance); % Optimal coefficients
        coef = [B0; B(:,idxLambdaMinDeviance)];


        hitLabel = [];
        
        Cvind =      crossvalind('Kfold',label,5);
        fprintf('\n');
        for r = 1:5
            fprintf('%d',r)
            trind = Cvind ~=r ;
            feat = Feature(trind,:);
            mdl =  fitglm(feat,label(trind),"Distribution","binomial",'Link','logit');

            [B,FitInfo] =lassoglm(feat,label(trind),"binomial", 'CV', 10);
        idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
        B0 = FitInfo.Intercept(idxLambdaMinDeviance); % Intercept
        coefficients = B(:, idxLambdaMinDeviance); % Optimal coefficients
        coef = [B0; B(:,idxLambdaMinDeviance)];


            teind = Cvind == r;
            [result,ci] = predict(mdl,Feature(teind,:));

            yhat = glmval(coef,Feature(teind,:),'logit');
            yhatBinom = (yhat>=0.5);

            acc(r) = length(find(result == label(teind)))/length(result);
        end
        fprintf('\n');
        CVAcc = mean(acc(r));
        fprintf('CV accuracy: %.2f\n',CVAcc);
        param.CVAcc = CVAcc;

        param.trD.mode = 'testing';
        
    catch
        fprintf('Training failed..!\n')
        hitLabel = -1;
    end
else
    try
        
        %         Feature = (Feature - repmat( param.trD.featparam(1), size(Feature,1),1))./repmat(param.trD.featparam(2), size(Feature,1),1);
          [~,sc] = predict(param.trD.mdl, Feature);
        if isfield(param,'DSP') || isfield(param.DTP,'W')
            if isfield(param.DSP,'W') || isfield(param.DTP,'W')
                
%                 for tr = 1:size(Feature,1)
%                     % check if ERP amplitude was 0 (ignored trial, only DSP w0 values) 
%                    if sum(Feature(tr,:) == reshape(param.DSP.w0(1:param.DSP.nf,:),1,size(Feature,2))) == size(Feature,2)
%                        sc(tr,:) = 0; 
%                    end
%                 end

                Nt = size(Feature,1)/(param.repeat*param.NumStims);
                tescore_bl = reshape(sc(:,end),param.repeat,param.NumStims,Nt);
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