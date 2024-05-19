function  param  = updateDSP(EP_1block,output_init,param)
lambda = param.DSP.lambda;
tarnarratio = 1/param.NumStims;
Sb_old = param.DSP.Sb;
Sw_old = param.DSP.Sw;
Swinv_old = param.DSP.Swinv;
mu_old = param.DSP.mu;
W_old = param.DSP.W;
N_old = param.DSP.N;


window = param.DSP.window;


erptar = EP_1block.nar(window,:,:,output_init);
erpnar = EP_1block.nar(window,:,:,setdiff(1:param.NumStims,output_init));
erpnar = erpnar(:,:,:);

% estimate lambda
mu_all = mu_old(:,:,1) - mu_old(:,:,2);
mu_all_new = mean(erptar,3,'omitnan')-mean(erpnar,3,'omitnan');
% lambda_tmp = [];
% for ch = 1:size(mu_all,2)
% lambda_tmp(ch) = norm(mu_all(:,ch) - mu_all_new(:,ch),2);
% end

for ch = 1:size(mu_all,2)
lambda_tmp(ch) = dot(mu_all(:,ch),mu_all_new(:,ch))/(norm(mu_all(:,ch))*norm(mu_all_new(:,ch)));
end
% % lambda = 1-exp(-0.05*mean(lambda_tmp));
% % lambda = mean(1-1./(1+lambda_tmp./20));
 
lambda = mean(0.1 * exp(-(lambda_tmp).^2 / (2 * 0.2^2)));

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
    S_i(:,:,s) = sum(S,3);%mean(S,3);
end
% Sw_new = (1-lambda)*Sw_old + lambda*sum(S_i,3,'omitnan');
Sw_new = Sw_old + sum(S_i,3,'omitnan');

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
param.DSP.lambda = lambda;
param.DSP.CosSim = lambda_tmp;
% param.DSP.Swinv = Swinv_new;
end




function invCov = updateInverseCovariance(invCov, x, lambda)
    % invCov: ?쁽?옱 怨듬텇?궛 ?뻾?젹?쓽 ?뿭?뻾?젹
    % x: ?깉濡? 異붽??맂 ?뜲?씠?꽣 ?룷?씤?듃 (?뿴 踰≫꽣)
    % lambda: ?뒪移쇰씪, 遺꾩궛 ?뾽?뜲?씠?듃?뿉 ?궗?슜?릺?뒗 ?뒪耳??씪 ?씤?옄

    % Sherman-Morrison-Woodbury 怨듭떇 ?쟻?슜
    u = x;
    v = x';
    invCov = (1/(1-lambda))*(invCov - (invCov * u * v * invCov) / ((1-lambda)/lambda + v * invCov * u));
end