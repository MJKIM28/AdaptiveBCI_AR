function [up,Posterior] = checkupdate2(score,pdfs,threshold,output)
% score: trial x stim


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

end
