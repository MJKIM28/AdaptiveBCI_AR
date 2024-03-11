function [Posterior] = getpostBayes(score,pdfs)
% score: trial x stim

%% Bayesian 
x = pdfs.x;
px_tar = pdfs.px_tar; % pdf of P(x|target);
px_nar = pdfs.px_nar; % pdf of P(x|nontarget);
Nstim = size(score,2);
Posterior = [];
for stim = 1:Nstim
    %-- mean score
    scoremean = mean(score(:,stim));
    %-- get x
    [~,xx] = min(abs(scoremean - x));
    pt = 1/Nstim; pn = (Nstim-1)/Nstim; % prior
    %-- update posterior
    px_t = px_tar(xx); % P(x|target)
    px_n = px_nar(xx); % P(x|nontarget)
    denom = (px_t*pt + px_n*pn);
    pt_x = (px_t*pt)/denom;
    pn_x = (px_n*pn)/denom;
    Posterior(stim,1) = pt_x; % posterior P(target|x)
    Posterior(stim,2) = pn_x; % posterior P(nontarget|x)
end
end