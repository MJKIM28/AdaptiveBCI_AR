close all; clear all

%%

LAMBDA1 = [0.01:0.01:0.1];
LAMBDA2 = [0.01:0.01:0.1];

%%

for lambda1 = LAMBDA1
    for lambda2 = LAMBDA2
SS03_main(lambda1,lambda2,0.3,0.5)
    end
end

% %%
% 
% Acc = squeeze(sum(outputs ==answers)./sum(~isnan(answers)))'
% mean(Acc)
% % figure;
% % for i = 1:Ntime; subplot(3,Ntime,i); topoplot(DSPs{1,1,s}(:,1,i),chanlocIn); colorbar; end
% % for i = 1:Ntime; subplot(3,Ntime,Ntime+i); topoplot(DSPs{end,1,s}(:,1,i),chanlocIn); colorbar;end
% % for i = 1:Ntime; subplot(3,Ntime,Ntime*2+i); topoplot(fixed.DSPs{s}(:,1,i),chanlocIn); colorbar;end
% [~,res1_] = max(outputs1);
% [~,res2_] = max(outputs2);
% 
% Acc1 = squeeze(sum(squeeze(res1_) ==answers)./sum(~isnan(answers)))'
% Acc2 = squeeze(sum(squeeze(res2_) ==answers)./sum(~isnan(answers)))'


