
function [C,param] = MWDSP(EP,param)

if strcmp(param.trD.mode,'training')
    winsize = 0.2;
    stepsize = 0.02;
    Win_temp = param.Baseline+1:param.Baseline+winsize*param.Fs;
    Win = bsxfun(@plus,Win_temp',0:stepsize*param.Fs:(param.Epocline-winsize*param.Fs));

    NtrAll = size(EP.tar,3);
    Nfit = fix(NtrAll/2);
    Nfit_ntar = Nfit*(param.NumStims-1);

    Nte = NtrAll-Nfit;
    Nte_nt = (NtrAll-Nfit)*(param.NumStims-1);

    Nch = size(EP.tar,2);
    STtar = EP.tar(:,:,1:Nfit);
    STntar = EP.nar(:,:,1:Nfit_ntar);

    STtar_te = EP.tar(:,:,Nfit+1:end);
    STntar_te = EP.nar(:,:,Nfit_ntar+1:end);
    param.MWDSP.nf  = 10;
    U1 = []; U2 = []; MeanX1 = []; MeanX2 = [];
    rho1 = []; rho2 = [];
    for ii = 1:size(Win,2)
        X{1} = STtar(Win(:,ii),:,:);
        X{2} = STntar(Win(:,ii),:,:);
        [W,w0,Sw,Sb,mu,N] = DSP(X,param.DSP.theta);

        param.MWDSP.DSP{ii}.W = W;
        param.MWDSP.DSP{ii}.w0 = w0;
        param.MWDSP.DSP{ii}.Sw = Sw;
        param.MWDSP.DSP{ii}.Sb = Sb;
        param.MWDSP.DSP{ii}.mu = mu;
        param.MWDSP.DSP{ii}.N = N;

        MeanX1(:,:,ii) = mean(X{1},3,'omitnan');
        MeanX2(:,:,ii) = mean(X{2},3,'omitnan');

        U1(:,:,ii) = param.MWDSP.DSP{ii}.W(:,1:param.MWDSP.nf)'*MeanX1(:,:,ii)' + param.MWDSP.DSP{ii}.w0(1:param.MWDSP.nf,:);
        U2(:,:,ii) = param.MWDSP.DSP{ii}.W(:,1:param.MWDSP.nf)'*MeanX2(:,:,ii)' + param.MWDSP.DSP{ii}.w0(1:param.MWDSP.nf,:);


        Xtar_te = permute(STtar_te(Win(:,ii),:,:),[2,1,3]); % CH x Time x Trial (Repeat X Block)
        Xtar_te = reshape(Xtar_te,Nch,length(Win(:,ii))*Nte);

        Xntar_te = permute(STntar_te(Win(:,ii),:,:),[2,1,3]); % CH x Time x Trial (Repeat X Block)
        Xntar_te = reshape(Xntar_te,Nch,length(Win(:,ii))*Nte_nt);


        Yt_test = param.MWDSP.DSP{ii}.W(:,1:param.MWDSP.nf)'*Xtar_te + repmat(param.MWDSP.DSP{ii}.w0(1:param.MWDSP.nf,:),1,Nte);
        Yn_test = param.MWDSP.DSP{ii}.W(:,1:param.MWDSP.nf)'*Xntar_te + repmat(param.MWDSP.DSP{ii}.w0(1:param.MWDSP.nf,:),1,Nte_nt);

        Yt = reshape(Yt_test,param.MWDSP.nf,length(Win(:,ii)),Nte);
        Yn = reshape(Yn_test,param.MWDSP.nf,length(Win(:,ii)),Nte_nt);


        UU1 = U1(:,:,ii);
        UU2 = U2(:,:,ii);
        for tt = 1:Nte
            Y = Yt(:,:,tt);

            rho1(tt,ii) = norm(UU1(:)-Y(:));
        end
        t_rho = Nte+1;
        for tt = 1:Nte_nt
            Y = Yn(:,:,tt);
            rho1(t_rho,ii) = norm(UU1(:)-Y(:));
            t_rho = t_rho + 1;
        end

        for tt = 1:Nte
            Y = Yt(:,:,tt);
            rho2(tt,ii) = norm(UU2(:)-Y(:));
        end
        t_rho = Nte+1;
        for tt = 1:Nte_nt
            Y = Yn(:,:,tt);
            rho2(t_rho,ii) = norm(UU2(:)-Y(:));
            t_rho = t_rho + 1;
        end
    end

    % yy = [ones(Nte,1);-ones(Nte_nt,1)];
    % [b,se,pval,finalmodel,stats] = stepwisefit(rho1-rho2,yy)

    yy = [ones(Nte,1);zeros(Nte_nt,1)];
    xx = rho1-rho2;

    mdl = stepwiseglm(xx, yy,'constant','upper','linear','distr','binomial');

    if mdl.NumEstimatedCoefficients > 1
    inmodel = [];
    for i=2:mdl.NumEstimatedCoefficients
        inmodel = [inmodel str2num(mdl.CoefficientNames{i}(2:end))];
    end
    else
        inmodel = 2:size(xx,2); 
    end

    param.MWDSP.TW.Win = Win;
    param.MWDSP.TW.sel = inmodel;

    param.MWDSP.TW.mdl = fitcsvm(xx(:,inmodel),yy);
    param.MWDSP.TW.mdl= incrementalLearner(param.MWDSP.TW.mdl);

    param.MWDSP.Template.U{1} = U1;
    param.MWDSP.Template.U{2} = U2;
    param.MWDSP.Template.Mean{1} = MeanX1;
    param.MWDSP.Template.Mean{2} = MeanX2;

    param.trD.mode = 'testing';
    C = [];
else

    ST = EP.nar;

    Nt = size(ST,5);
    rho_i = 1;
    rho1 = []; rho2 = [];
    for ii = param.MWDSP.TW.sel
        Temp = permute(ST(param.MWDSP.TW.Win(:,ii),:,:,:,:),[2,1,3,4,5]); % CH x Time x Repeat x Stim X Block
        Temp = reshape(Temp,size(ST,2),length(param.MWDSP.TW.Win(:,ii))*param.repeat*param.NumStims*Nt);
        Yt_test = param.MWDSP.DSP{ii}.W(:,1:param.MWDSP.nf)'*Temp...
            + repmat(param.MWDSP.DSP{ii}.w0(1:param.MWDSP.nf,:),1,param.repeat*param.NumStims*Nt);

        V = reshape(Yt_test,param.MWDSP.nf,length(param.MWDSP.TW.Win(:,ii)),param.repeat*param.NumStims*Nt);

        UU1 = param.MWDSP.Template.U{1}(:,:,ii);
        UU2 = param.MWDSP.Template.U{2}(:,:,ii);
        for tt = 1:size(V,3)
            Y = V(:,:,tt);
            rho1(tt,rho_i) = norm(UU1(:)-Y(:));
        end
        for tt = 1:size(V,3)
            Y = V(:,:,tt);
            rho2(tt,rho_i) = norm(UU2(:)-Y(:));
        end
        rho_i = rho_i+1;

    end
    xx = rho1-rho2;
    [k,sc] = predict(param.MWDSP.TW.mdl,xx);

    [~,C] = max(sum(reshape(sc(:,end),param.repeat,param.NumStims)));
end