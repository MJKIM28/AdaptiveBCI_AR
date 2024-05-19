function mdlnew = getLDA(Feat,label)
% label 0: nontarget
% label 1: target

classtype = unique(label);

K = length(classtype);
[N,D] = size(Feat);

mdlnew.Mu = zeros(K,size(Feat,2));
X = cell(K,1);
sigma = zeros(D,D,K);
for k = 1:K
    ind = find(label==classtype(k));
    X{k} = Feat(ind,:);
    M = mean(X{k},'omitnan');

    
    mdlnew.N(k) = length(ind);
    mdlnew.Mu(k,:) = M;

    sigma(:,:,k) = (X{k} - mdlnew.Mu(k,:))'*(X{k} - mdlnew.Mu(k,:));
end
mdlnew.MuTotal = mean(mdlnew.Mu,1);
mdlnew.Sigma = sum(sigma,3)/(N-K);
mdlnew.invSigma = inv(mdlnew.Sigma);

mdlnew.Coeffs = getLDAcoeffs(mdlnew);

