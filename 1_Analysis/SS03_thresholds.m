Thresholds = [0.1:0.2:1];
SubNames =  {'Subtest12','Subtest13','Subtest14','Subtest15','Subtest17',...
    'Subtest18','Subtest19','Subtest20','Subtest21','Subtest22'};
CondName = {'Pre','Main 1','Main 2','Main 3','Main 4','Main 5','Post'};
%%

AccAll = NaN([size(Thresholds),7,10]);
for l1 = 1:length(Thresholds)
    for l2 = 1:length(Thresholds)
filename = ['adaptiveSPHDCA\SS03_exp_adaptiveSPHDCA_LDA_l1_0.03_l2_0.03_thres1_',num2str(Thresholds(l1)),'_thres2_',num2str(Thresholds(l2)),'.mat'];

load(filename)

Acc = squeeze(sum(outputs ==answers)./sum(~isnan(answers)))';
AccAll(l1,l2,:,:) = Acc';

fprintf('*')
    end
    fprintf('.')
end
fprintf('\n')

save('adaptiveSPHDCA\Acc_thresholds','AccAll')
%%
load('adaptiveSPHDCA\Acc_thresholds','AccAll')


figure;
for i = 1:size(AccAll,4) % subject
for j = 1:size(AccAll,3) % session
    subplot(10,size(AccAll,3),(i-1)*7+j)
    imagesc(AccAll(:,:,j,i))
    caxis([0.6 1])
    set(gca,'xtick',1:length(Thresholds),'xticklabel',Thresholds,'ytick',1:length(Thresholds), ...
        'yticklabel',Thresholds)
    colormap jet
if i == 1
    title(CondName{j})
end
end
end

figure;
for j = 1:size(AccAll,3)
    subplot(1,size(AccAll,3),j)
    imagesc(mean(AccAll(:,:,j,:),4,'omitnan'))
%     caxis([0.6 1])
colorbar
    set(gca,'xtick',1:length(Thresholds),'xticklabel',Thresholds,'ytick',1:length(Thresholds), ...
        'yticklabel',Thresholds)
    colormap jet
if i == 1
    title(CondName{j})
end
end

figure;
for i = 1:size(AccAll,4)
    subplot(1,10,i)
    imagesc(mean(AccAll(:,:,:,i),3,'omitnan'));
    colorbar
    set(gca,'xtick',1:length(Thresholds),'xticklabel',Thresholds,'ytick',1:length(Thresholds), ...
        'yticklabel',Thresholds)
    colormap jet
end

figure;
imagesc(mean(mean(AccAll,3,'omitnan'),4,'omitnan'));
   set(gca,'xtick',1:length(Thresholds),'xticklabel',Thresholds,'ytick',1:length(Thresholds), ...
        'yticklabel',Thresholds)
%%
start_color = [0, 0, 255]; % Blue
end_color = [255, 0, 0]; % Red

% Number of samples
n_samples = 5;

% Create the gradient
gradient = zeros(n_samples, 3); % Initialize the matrix for RGB values
for i = 1:3
    gradient(:, i) = linspace(start_color(i), end_color(i), n_samples);
end

figure; 
for th1 = 1:5
    subplot(5,1,th1);hold on; 
for th2 = 1:5
plot(squeeze(mean(AccAll(th1,th2,:,6:end),4,'omitnan'))','color',gradient(th2,:)/255)
end
end
