function mdlnew = updateSingleClassfier(Feat,mdl,lambda)
% Pooled mean

% Sigma = mdl.Sigma;
% invSigma = mdl.invSigma;
% N = mdl.N;

Meanold = mdl.MuTotal;
MeanNow = mean(Feat);
lambda_tmp =  dot(Meanold,MeanNow)/(norm(Meanold)*norm(MeanNow));

% % lambda = mean(1-1./(1+lambda_tmp./20));
 lambda = mean(0.06 * exp(-(lambda_tmp).^2 / (2 * 0.2^2)));

%%
Meannew = (1-lambda)*Meanold+ lambda*MeanNow;
% covEst = ((n - 2) / (n - 1)) * Sigma + (X - oldMean) * (X - meanEst)' / n;
% invSigmanew = updateInverseCovariance(mdl.invSigma, Meannew', lambda);
invSigmanew = mdl.invSigma;
%%
mdlnew = mdl;
mdlnew.invSigma = invSigmanew;
Diff = -diff(mdlnew.Mu,1);
w = mdlnew.invSigma*Diff';
mdlnew.Coeffs(1).Linear = w;
mdlnew.Coeffs(2).Linear = -w;

mdlnew.MuTotal = Meannew;

mdlnew.Coeffs(1).Const = -mdlnew.Coeffs(1).Linear'*mdlnew.MuTotal'+log(mdl.N(1)/mdl.N(2));
mdlnew.Coeffs(2).Const = -mdlnew.Coeffs(1).Const;
mdlnew.lambda = lambda;
mdlnew.CosSim = lambda_tmp;
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