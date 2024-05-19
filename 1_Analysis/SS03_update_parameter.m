
lambdas =  0.01:0.01:0.1;
SubNames = {'Subtest17','Subtest18','Subtest19','Subtest20','Subtest21','Subtest22'};
CondName = {'Pre','Main 1','Main 2','Main 3','Main 4','Main 5','Post'};

AccAll = NaN([size(lambdas),7,6]);
for l1 = 1:length(lambdas)
    for l2 = 1:length(lambdas)
filename = ['adaptiveSPHDCA\SS03_5_MWDCPM_cond_comb3_adaptiveSPHDCA_LDA_l1_',num2str(lambdas(l1)),'_l2_',num2str(lambdas(l2)),'.mat'];

load(filename)

Acc = squeeze(sum(outputs ==answers)./sum(~isnan(answers)))';
AccAll(l1,l2,:,:) = Acc';

fprintf('*')

    end
        fprintf('.')

end

save('adaptiveSPHDCA\Acc_updateparam','AccAll')

%%
load('adaptiveSPHDCA\Acc_updateparam')
figure;
for i = 1:size(AccAll,3)
    subplot(1,size(AccAll,3),i)
    imagesc(mean(AccAll(:,:,i,:),4,'omitnan'))
%     caxis([0.5 1])
    set(gca,'xtick',1:length(lambdas),'xticklabel',lambdas,'ytick',1:length(lambdas), ...
        'yticklabel',lambdas)
    colorbar

end


figure;
imagesc(mean(mean(AccAll,3,'omitnan'),4,'omitnan'));
   set(gca,'xtick',1:length(lambdas),'xticklabel',lambdas,'ytick',1:length(lambdas), ...
        'yticklabel',lambdas)

figure;
for i = 1:size(AccAll,4)
for j = 1:size(AccAll,3)
    subplot(6,size(AccAll,3),(i-1)*7+j)
    imagesc(AccAll(:,:,j,i))
    caxis([0.6 1])
    set(gca,'xtick',1:length(lambdas),'xticklabel',lambdas,'ytick',1:length(lambdas), ...
        'yticklabel',lambdas)
    colormap jet
if i == 1
    title(CondName{j})
end
end

end
colorbar
xlabel('LDA')
ylabel('DSP')
