function [up_check,upids,Posterior,trs,Dist] = checkupdate(score,pdfs,threshold,output,type)
% score: trial x stim
[Ntrial,Nstim] = size(score);

%% Bayesian 
% x = pdfs.x;
% px_tar = pdfs.px_tar; % pdf of P(x|target);
% px_nar = pdfs.px_nar; % pdf of P(x|nontarget);
% Nstim = size(score,2);
% Posterior = [];
% for stim = 1:Nstim
%     %-- mean score
%     scoremean = mean(score(:,stim));
%     %-- get x
%     [~,xx] = min(abs(scoremean - x));
%     pt = 1/Nstim; pn = (Nstim-1)/Nstim; % prior
%     %-- update posterior
%     px_t = px_tar(xx); % P(x|target)
%     px_n = px_nar(xx); % P(x|nontarget)
%     denom = (px_t*pt + px_n*pn);
%     pt_x = (px_t*pt)/denom;
%     pn_x = (px_n*pn)/denom;
%     Posterior(stim,1) = pt_x; % posterior P(target|x)
%     Posterior(stim,2) = pn_x; % posterior P(nontarget|x)
% end
% 
% if max(Posterior(:,1)) > threshold
%     up = true;
% else
%     up = false;
% end

%% Mahalanobis distance
% if strcmp(type,'md')
Dist = [];
for ii = 1:Nstim
    Dist(ii) = mahal(mean(score(:,ii)),reshape(score(:,setdiff(1:Nstim,ii)),Ntrial*(Nstim-1),1));
end
[maxval,maxid] = max(Dist);
if maxid == output %&& maxval> threshold 
    up_check = true;
else
    up_check = false;
end
% Posterior = maxval;

%% Z scored score
% MM = mean(score(:));
% STD = std(score(:));
% % ZZ = (mean(score) - MM)/STD;
% ZZ = (score - MM)/STD; % for averaged feature
% [maxval,maxid] = max(ZZ);
% if maxval> threshold && maxid == output
%     up = true;
% else
%     up = false;
% end
% Posterior = maxval;

%% ratio

% % elseif strcmp(type,'ratio')
% 
% [Sorted,idx] = sort(score,2,"descend");
% 
% trs = find(idx(:,1) == output); % max score == prediction
% Posterior = 1 - Sorted(:,2)./Sorted(:,1);
% 
% ids = Posterior(trs) > threshold;
% upids = trs(ids);
% % end
upids = 0; Posterior = 0; trs = 0;
end
