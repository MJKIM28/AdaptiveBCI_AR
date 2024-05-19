%% recognition test

%%
fpath = 'E:\1_연구\분석코드\2024_adaptive_BCI\Data\Pilot';

SubName = {'Subtest17','Subtest18','Subtest19','Subtest20','Subtest21','Subtest22'};
Nsub = length(SubName);

Nsessmain = 5;
%%
ANSWER = [ones(5,1);zeros(5,1)];
Acc_rec = NaN(Nsub,Nsessmain);
for s = 1:Nsub
    load([fpath,'\',SubName{s}])

    for sess = 1:Nsessmain

        if ~isempty(vars.Recog{sess})

            for i = 1:length(vars.Recog{sess})
                if isempty(vars.Recog{sess}(i).response)
                    vars.Recog{sess}(i).response = NaN;
                end

               
            end

             Resp = [vars.Recog{sess}.response];
                Img = [vars.Recog{sess}.image];
                [~,id] = sort(Img,'ascend');
                Acc_rec(s,sess) = mean(Resp(id)' == ANSWER,'omitnan');

        end
    end
end
Acc_rec = Acc_rec*100;
save('Accuracy_recog.mat','Acc_rec')
%%
load('Accuracy_recog.mat')
Names = {'Pilot16','Pilot17','Pilot18','Pilot19','Pilot20','Pilot21'};


% % Define the start and end colors
% start_color = [0, 0, 255]; % Blue
% end_color = [255, 0, 0]; % Red
% 
% gradient = zeros(10, 3); % Initialize the matrix for RGB values
% for i = 1:3
%     gradient(:, i) = linspace(start_color(i), end_color(i), 10);
% end
colors = [0 0 255;255 0 0;0 128 0;128 0 128;0 128 128;169 169 169]./255;

figure;
a = [];
plot(mean(Acc_rec,'omitnan'),'Color','black','marker','s','linewidth',3)
for s = 1:Nsub
    hold on;
    a(s) = plot(Acc_rec(s,:)'*100,'Color',colors(s,:),'marker','o','LineWidth',2,'MarkerFaceColor',colors(s,:));    

xlim([0.5 5.5])
ylim([40 105])
set(gca,'XTick',1:7,'XTickLabel',{'Main 1','Main 2','Main 3','Main 4','Main 5'},'FontSize',15)
ylabel('Accuracy (%)')
end
legend(a,Names,'Location','eastoutside')
% set(gcf,'Position',[ 329   283   1318   697])


%% 
load('Accuracy_BCI_exp.mat')


figure;
for s = 1:Nsub
    subplot(2,3,s);
    hold on;
    a1 = plot(Acc_fix(s+4,2:end-1)','Color',[0.6 0.6 0.6],'marker','o','LineWidth',2,'MarkerFaceColor',[0.6 0.6 0.6]);    
    
    a0 = plot(Acc_ad_selbl(s+4,2:end-1)','Color','k','marker','o','LineWidth',2,'MarkerFaceColor','k');    
    
    a3 = plot(Acc_rec(s,:)','Color','b','marker','diamond','LineWidth',2,'MarkerFaceColor','b');
    

xlim([0.5 5.5])
ylim([0 105])
set(gca,'XTick',1:7,'XTickLabel',{'Main 1','Main 2','Main 3','Main 4','Main 5'},'FontSize',15)
ylabel('Accuracy (%)')
title(Names{s})
end
legend([a1 a0 a3],{'BCI (fix)','BCI (Ad)','Recognition test'},'Location','eastoutside')
% set(gcf,'Position',[ 329   283   1318   697])

ac = Acc_ad_selbl(5:end,2:end-1);
[r,p] = corr(ac(2,:)',Acc_rec(2,:)',"rows","pairwise")