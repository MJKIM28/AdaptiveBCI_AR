function [W, w0,Sw_reg,Sb, mu,N] = DSP(X,theta)
% X{i}: class i, time x channel x trial
% theta: the regularization level of Sw
% Assume that the time length and the number of channel is same for every classes

Nclass = length(X);
[T,C,~] = size(X{1});

S = zeros(C,C,Nclass);
mu = zeros(T,C,Nclass);
for i = 1:Nclass
    N(i) = size(X{i},3);
    
    mu(:,:,i) = mean(X{i},3);
    
    for t = 1:N(i)
        Wd = X{i}(:,:,t) - mu(:,:,i);
        S_i = Wd'*Wd;
        S(:,:,i) = S(:,:,i) + S_i;
    end
    S(:,:,i) = S(:,:,i)/N(i);
end

Sw = sum(S,3);
Sw_reg = (1-theta)*Sw + theta*eye(size(Sw)); % regularize Sw with regularization level theta

mu_all = mean(mu.*permute(repmat(N/sum(N),T,1,C),[1,3,2]),3);

Sb = 0;
for i = 1:Nclass
    Bd = mu(:,:,i) - mu_all;
    Sb = Sb + N(i)*(Bd'*Bd);
end

[W, U] = eig(inv(Sw_reg)*Sb);

w0 = -W'*mu_all';

end