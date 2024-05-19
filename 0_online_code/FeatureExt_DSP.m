function [Feature,labels, param] = FeatureExt_DSP(EP, param)
% Nt = size(EP.lat,3);
window = param.DSP.window;

if strcmp(param.trD.mode,'training')
    Nt = size(EP.tar,3);
    Nntr = size(EP.nar,3);

    STtar = EP.tar;
    STntar = EP.nar;
    
    X{1} = STtar(window,:,:);
    X{2} = STntar(window,:,:);
    
    %-- updated: to get Sw,Sb,mu (@20221025)
    if ~isfield(param.DSP,'W') 
        [W,w0,Sw,Sb,mu,N,Swinv] = DSP(X,param.DSP.theta);
        param.DSP.W = W;
        param.DSP.w0 = w0;
        param.DSP.Sw = Sw;
        param.DSP.Sb = Sb;
        param.DSP.mu = mu;
        param.DSP.N = N;
        param.DSP.Swinv = Swinv;
        
    else
        fprintf('use previous DSP weight\n')
    end
    
    Yt = [];
    Yn = [];
    for t =1:Nt
        Yt(:,:,t) = param.DSP.W'*STtar(window,:,t)' + param.DSP.w0;
    end
    for t = 1:Nntr
        Yn(:,:,t) = param.DSP.W'*STntar(window,:,t)' + param.DSP.w0;
    end
    
    %% get fisher score and select feature
    if strcmp(param.DSP.type,'fscore')
        sizetar = size(Yt,3);
        sizenar = size(Yn,3);
        fscore = [];
        for c = 1:size(Yt,1)
            for ss = 1:length(window)
                MeanAll = mean(cat(3,Yt(c,ss,:),Yn(c,ss,:)),3);
                nom = sizetar*(mean(Yt(c,ss,:),3) - MeanAll).^2 + sizenar*(mean(Yn(c,ss,:),3) - MeanAll).^2;
                denom = sizetar*std(Yt(c,ss,:),[],3).^2 + sizenar*std(Yn(c,ss,:),[],3).^2;
                fscore(c,ss) = nom/denom;
            end
        end
        
        featsel = fscore > mean(fscore(:));
        Compsel = find(sum(featsel,2) > length(window)*0.1);
        featsel(setdiff(1:size(Yt,1),Compsel),:) = 0;
        featid = find(featsel);
        
        traintar = []; trainnar = [];
        for t = 1:sizetar
            traintmp1 = Yt(:,:,t);
            traintar(t,:) = traintmp1(featid);
        end
        for tt = 1:sizenar
            traintmp2 = Yn(:,:,tt);
            trainnar(tt,:) = traintmp2(featid);
        end
        param.DSP.Featid = featid;
        
        Feature = [traintar; trainnar];
        labels = [ones(sizetar,1); -ones(sizenar,1)];
        
        
    elseif strcmp(param.DSP.type,'nf')
        
        traintmp1 = reshape(Yt(1:param.DSP.nf,:,:),param.DSP.nf*size(Yt,2),size(Yt,3));
        traintmp2 = reshape(Yn(1:param.DSP.nf,:,:),param.DSP.nf*size(Yn,2),size(Yn,3));
        %
        Feature = [traintmp1 traintmp2]';
        labels = [ones(size(traintmp1,2),1); -ones(size(traintmp2,2),1)];
    end
else
    
    
    
    %      Yt_test = [];
    %     for t = 1:Nt*param.NumStims*param.repeat
    %         Yt_test(:,:,t) = param.DSP.W'*ST(window,:,t)' + param.DSP.w0;
    %     end
    ST = EP.nar;
    
    Nt = size(ST,5);
    
    Temp = permute(ST(window,:,:,:,:),[2,1,3,4,5]); % CH x Time x Repeat x Stim X Block
    Temp = reshape(Temp,size(ST,2),length(window)*param.repeat*param.NumStims*Nt);
    Yt_test = param.DSP.W'*Temp + repmat(param.DSP.w0,1,param.repeat*param.NumStims*Nt);
    
    %% use selected feature
    if strcmp(param.DSP.type,'fscore')
        testdata = [];
        testtmp = reshape(Yt_test,size(Yt_test,1),length(window),param.repeat*param.NumStims*Nt);
        for ttt = 1:size(testtmp,3)
            tmp = testtmp(:,:,ttt);
            testdata(ttt,:) = tmp(param.DSP.Featid);
        end
        Feature = testdata;
    
    elseif strcmp(param.DSP.type,'nf')
        testtmp = reshape(Yt_test(1:param.DSP.nf,:,:),param.DSP.nf*length(window),param.repeat*param.NumStims*Nt);
        Feature = testtmp';
    end
    
    labels = [];
    
    
end