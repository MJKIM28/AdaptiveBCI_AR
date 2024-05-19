function prob = pred_mj(mdl,X)
K = 2;
[N,D] = size(X);
probabilities = zeros(1, K);

% Sigma = mdl.Sigma; % 공분산 행렬
% invSigma = mdl.invSigma;
% for k = 1:K
%     mu_k = mdl.Mu(k,:);
%     Diff = X - mu_k;
%     likelihood = exp(-0.5 * (Diff*invSigma) * Diff') / (sqrt((2*pi)^(d) * det(Sigma)));
%     probabilities(k) = mdl.Prior(k) * likelihood;
% end
% prob = probabilities / sum(probabilities);



scores = zeros(N,K);
for k = 1:K
    linearCoeff = mdl.Coeffs(1, k).Linear;
    constCoeff = mdl.Coeffs(1, k).Const;
    scores(:,k) = X*linearCoeff + constCoeff;
end


prob = 1./(1+exp(-scores));
