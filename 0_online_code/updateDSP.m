function  param  = updateDSP(EP_1block,output_init,param)
lambda = param.DSP.lambda;
tarnarratio = 1/param.NumStims;
Sb_old = param.DSP.Sb;
Sw_old = param.DSP.Sw;
mu_old = param.DSP.mu;
W_old = param.DSP.W;
N_old = param.DSP.N;


window = param.Baseline + 1:param.winsize:param.Baseline+param.Epocline;


erptar = EP_1block.nar(window,:,:,output_init);
erpnar = EP_1block.nar(window,:,:,setdiff(1:param.NumStims,output_init));
erpnar = erpnar(:,:,:);

%-- upated mu
mu_new(:,:,1) = (1-lambda)*mu_old(:,:,1) + lambda*mean(erptar,3,'omitnan');
mu_new(:,:,2) = (1-lambda)*mu_old(:,:,2) + lambda*mean(erpnar,3,'omitnan');

%-- update Sw
Wd{1} = erptar - mu_new(:,:,1);
Wd{2} = erpnar - mu_new(:,:,2);
for s = 1:length(Wd)
    N(s) = size(Wd{s},3);
    N_new(s) = N_old(s)+N(s);
    for i = 1:N(s)
        S(:,:,i) = Wd{s}(:,:,i)'*Wd{s}(:,:,i);
    end
    S_i(:,:,s) = mean(S,3);
end
Sw_new = (1-lambda)*Sw_old + lambda*sum(S_i,3,'omitnan');

%-- update Sw^(-1)
Sw_inv = inv(Sw_new+0.00001*eye(size(Sw_new)));

%-- update Sb
Bd(:,:,1) = (1-tarnarratio)*(mu_new(:,:,1) - mu_new(:,:,2));
Bd(:,:,2) = (tarnarratio)*(mu_new(:,:,2) - mu_new(:,:,1));
Sb_new = 0;
for s = 1:2
    Sb_new =  Sb_new + Bd(:,:,s)'*Bd(:,:,s)*(N_new(s));
end

%-- update W
Sb_n = Sb_new;
W_new = [];
for n = 1:size(erptar,2)%param.DSP.nf
    if n > 1
        %--
        nume = Sb_n*W_new(:,n-1)*W_new(:,n-1)';
        denom = W_new(:,n-1)'*Sb_n*W_new(:,n-1); 
        Sb_n =(eye(size(Sb_n)) - nume/denom)* Sb_n;
    end
    num = W_old(:,n)'*Sw_new*W_old(:,n);   % Sw in numerator (fixed in 20240308)
    den = W_old(:,n)'*Sb_n*W_old(:,n); % Sb in denominator
    W_new(:,n) = (num/den)*Sw_inv*Sb_n*W_old(:,n);
    W_new(:,n) = W_new(:,n)/norm(W_new(:,n));

end
% [W_new, U] = eig(Sw_inv*Sb_new);

w0 = -W_new'*(tarnarratio*mu_new(:,:,1)+(1-tarnarratio)*mu_new(:,:,2))';

fprintf('DSP weight updated\n')
param.DSP.W = W_new;
param.DSP.Sw = Sw_new;
param.DSP.Sb = Sb_new;
param.DSP.mu = mu_new;
param.DSP.N = N_new;
param.DSP.w0 = w0;
end