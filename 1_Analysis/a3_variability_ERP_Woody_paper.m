close all; clear all
%% parameters
fpath = 'E:\1_연구\분석코드\2024_adaptive_BCI\Data\Pilot';

SubNameList1 = {'Subtest02','Subtest03','Subtest04','Subtest05'};

SubNameList2 = {'Subtest12','Subtest13','Subtest14','Subtest15','Subtest17',...
    'Subtest18','Subtest19','Subtest20','Subtest21','Subtest22'};

SubName = SubNameList2;

Nsub = length(SubName);

Fs = 500;
Nch = 31;
Tbase = 0.2;
Tepoc = 0.6;
Lepoc = (Tbase+Tepoc)*Fs;

Ntr_tr = 15;
Ntr_con  = 15;
Nsess_main = 5;
Nsess = Nsess_main  + 2;
Nrepeat = 10;
chanloc = readlocs('D:\[1]EEGBCI\[2]Research\Code&Algorithm\BP_channelLocs.locs');

ChList = {'Fp1'    'Fpz'    'Fp2'    'F7'    'F3'    'Fz'    'F4'    'F8' ...
    'FT9'    'FC5'    'FC1'    'FC2'    'FC6'    'FT10'    'T7'    'C3'...
    'Cz'    'C4'    'T8'    'CP5'    'CP1'    'CP2'    'CP6'    'P7'...
    'P3'    'Pz'    'P4'    'P8'    'O1'    'Oz'    'O2'};
Times = -0.2:1/Fs:0.6-1/Fs;
Ntime = length(Times);
CondName2 = {'Control','1st','2nd','3rd','4th','5th'};

load("TrialUsed.mat")
load("TrialUsed2.mat")
trial_use = [Trials(9:end); Trials2];


%% per component, use training mean as template
load(['ERP\value_ERPComp_FC1FC2.mat'])

Latencies = cat(1,P2lat,P3lat,N2lat);

CompName= {'P3a','P3b','N2'};
chdrawName{1} = {'Fz','FC1','FC2','Cz'};
chdrawName{2} = {'P3','Pz','P4'};
chdrawName{3} = {'O1','Oz','O2'};

Ncomp = size(chdrawName,2);

load(['ERP\ERP_tar'])

%% Woody's (initial template from training ERP, use separate time window)

Win{1} = find(Times >= 0.1 & Times <= 0.3);
Win{2} = find(Times >= 0.25 & Times <= 0.45);

% TimeWindow{1} = Win{1}; % P3a
% TimeWindow{2} = Win{2}; % P3b
% TimeWindow{3} = Win{1}; % N2

TimeWindow{1} = find(Times > 0 & Times < 0.6); % P3a
TimeWindow{2} = find(Times > 0 & Times < 0.6);  % P3b
TimeWindow{3} = find(Times > 0 & Times < 0.6);  % N2

limits = 50; % 100ms
threshold = 1; % 0.0001: iteration, 1:not iterate
LAGs = NaN(Ntr_con*Nrepeat,Ncomp,Nsess,Nsub);
STDs = [];
ERPalmean = cell(3,1);
ERPormean = cell(3,1);
for comp = 1:Ncomp
    ERPalmean{comp} = NaN(length(Times),length(chdrawName{comp}),Nsess,Nsub);
    ERPormean{comp} = NaN(length(Times),length(chdrawName{comp}),Nsess,Nsub);
end
ERPalsingle = cell(3,1);
ERPorsingle = cell(3,1);
for comp = 1:Ncomp
    ERPalsingle{comp} = NaN(length(Times),Ntr_con*Nrepeat,length(chdrawName{comp}),Nsess,Nsub);
    ERPorsingle{comp} = NaN(length(Times),Ntr_con*Nrepeat,length(chdrawName{comp}),Nsess,Nsub);
end
%%
for s = 1:Nsub
    fprintf('%d\n',s)

    chlist = find(~isnan(EP_tar_main(1,:,1,2,s)));

    for comp = 1:Ncomp
        chdraw = find(ismember(ChList,chdrawName{comp}));

        chids = find(ismember(chdraw,chlist));
        chdraw_use = chdraw(chids);
        chdraw_indv =  find(ismember(chlist,chdraw_use));

        e_tr = load([fpath,'\Train\Epoch\',SubName{s}]);

        %-- target
        Template = squeeze(mean(e_tr.Epoch.tar(TimeWindow{comp},chdraw_indv,:),3,'omitnan'));
        %% pre
trial_use = find(~isnan(EP_tar_pre(1,1,:,s)));

        if TimeWindow{comp}(1) - limits*Fs/1000 > 0 % signal after the window
            ERPouttemp1 = EP_tar_pre(TimeWindow{comp}(1) - limits*Fs/1000: TimeWindow{comp}(1)-1,chdraw_use,trial_use,s);
        else % signal from the first point
            ERPouttemp1 = EP_tar_pre(1: limits*Fs/1000,chdraw_use,trial_use,s);
        end
        if TimeWindow{comp}(end) + limits*Fs/1000 <= Ntime % signal before the window
            ERPouttemp2 = EP_tar_pre(TimeWindow{comp}(end) + 1: TimeWindow{comp}(end) + limits*Fs/1000,chdraw_use,trial_use,s);
        else % signal to the end point
            ERPouttemp2 = EP_tar_pre(TimeWindow{comp}(end) - limits*Fs/1000 + 1: end,chdraw_use,trial_use,s);
        end
        ERPout{1} = ERPouttemp1;
        ERPout{2} = ERPouttemp2;

        ERPs = EP_tar_pre(TimeWindow{comp},chdraw_use,trial_use,s);
        [ERPaligned, Lags,ERPtrhistory,Lagtrhistory,CORs] = Woodys(ERPs,limits,threshold,Fs,Template,ERPout);

        STDs(comp,1,s) = std(Lags(:)/Fs*1000);
        LAGs(trial_use,comp,1,s) = Lags(:);


        ERPsLong = EP_tar_pre(:,chdraw_use,:,s);
        ERPalignedLong = [];
        for tt = 1:length(trial_use)
            L = LAGs(tt,comp,1,s);
            if L < 0
                ERPalignedLong(:,:,tt) = cat(1,ERPsLong(abs(L)+1:end,:,tt), ERPsLong(1:abs(L),:,tt));
            elseif L>=0
                ERPalignedLong(:,:,tt) = cat(1,ERPsLong(end-abs(L)+1:end,:,tt), ERPsLong(1:end-abs(L),:,tt));
            end

        end
        ERPalsingle{comp}(:,trial_use,ismember(chdraw,chdraw_use),1,s) = permute(ERPalignedLong,[1 3 2]);
        ERPorsingle{comp}(:,trial_use,ismember(chdraw,chdraw_use),1,s) = permute(ERPsLong,[1 3 2]);
        ERPalmean{comp}(:,ismember(chdraw,chdraw_use),1,s) = mean(ERPalignedLong,3);
        ERPormean{comp}(:,ismember(chdraw,chdraw_use),1,s) = mean(ERPsLong,3);

        %% main
        for sess = 1:Nsess_main
            
            trial_use = find(~isnan(EP_tar_main(1,1,:,sess,s)));
            if s> 4 && trial_use(end) == 150
                trial_use(find(trial_use == 141):end) = [];
            end
            if ~isempty(trial_use)

            if TimeWindow{comp}(1) - limits*Fs/1000 > 0 % signal after the window
                ERPouttemp1 = EP_tar_main(TimeWindow{comp}(1) - limits*Fs/1000: TimeWindow{comp}(1)-1,chdraw_use,trial_use,sess,s);
            else % signal from the first point
                ERPouttemp1 = EP_tar_main(1: limits*Fs/1000,chdraw_use,:,sess,s);
            end
            if TimeWindow{comp}(end) + limits*Fs/1000 <= Ntime % signal before the window
                ERPouttemp2 = EP_tar_main(TimeWindow{comp}(end) + 1: TimeWindow{comp}(end) + limits*Fs/1000,chdraw_use,trial_use,sess,s);
            else % signal to the end point
                ERPouttemp2 = EP_tar_main(TimeWindow{comp}(end) - limits*Fs/1000 + 1: end,chdraw_use,trial_use,sess,s);
            end
            ERPout{1} = ERPouttemp1;
            ERPout{2} = ERPouttemp2;

            %%
            ERPs = EP_tar_main(TimeWindow{comp},chdraw_use,trial_use,sess,s);
            [ERPaligned, Lags,ERPtrhistory,Lagtrhistory,CORs] = Woodys(ERPs,limits,threshold,Fs,Template,ERPout);
            %%
            %                         figure; subplot(121); plot(mean(ERPs,3));
            %                         subplot(122); plot(mean(ERPaligned,3));

            STDs(comp,sess+1,s) = std(Lags(:)/Fs*1000);
            LAGs(trial_use,comp,sess+1,s) = Lags(:);

            %%
            ERPsLong = EP_tar_main(:,chdraw_use,trial_use,sess,s);
            ERPalignedLong = [];
            for tt = 1:length(trial_use)
                L = LAGs(tt,comp,sess+1,s);
                if L < 0
                    ERPalignedLong(:,:,tt) = cat(1,ERPsLong(abs(L)+1:end,:,tt), ERPsLong(1:abs(L),:,tt));
                elseif L>=0
                    ERPalignedLong(:,:,tt) = cat(1,ERPsLong(end-abs(L)+1:end,:,tt), ERPsLong(1:end-abs(L),:,tt));
                end

            end
            ERPalsingle{comp}(:,trial_use,ismember(chdraw,chdraw_use),sess+1,s) = permute(ERPalignedLong,[1 3 2]);
            ERPorsingle{comp}(:,trial_use,ismember(chdraw,chdraw_use),sess+1,s) = permute(ERPsLong,[1 3 2]);
            ERPalmean{comp}(:,ismember(chdraw,chdraw_use),sess+1,s) = mean(ERPalignedLong,3);
            ERPormean{comp}(:,ismember(chdraw,chdraw_use),sess+1,s) = mean(ERPsLong,3);
            end
        end
            %% post
            
trial_use = find(~isnan(EP_tar_post(1,1,:,s)));
if s> 4 && trial_use(end) == 150
                trial_use(find(trial_use == 141):end) = [];
            end
            if ~isempty(trial_use)

            if TimeWindow{comp}(1) - limits*Fs/1000 > 0 % signal after the window
                ERPouttemp1 = EP_tar_post(TimeWindow{comp}(1) - limits*Fs/1000: TimeWindow{comp}(1)-1,chdraw_use,trial_use,s);
            else % signal from the first point
                ERPouttemp1 = EP_tar_post(1: limits*Fs/1000,chdraw_use,trial_use,s);
            end
            if TimeWindow{comp}(end) + limits*Fs/1000 <= Ntime % signal before the window
                ERPouttemp2 = EP_tar_post(TimeWindow{comp}(end) + 1: TimeWindow{comp}(end) + limits*Fs/1000,chdraw_use,trial_use,s);
            else % signal to the end point
                ERPouttemp2 = EP_tar_post(TimeWindow{comp}(end) - limits*Fs/1000 + 1: end,chdraw_use,trial_use,s);
            end
            ERPout{1} = ERPouttemp1;
            ERPout{2} = ERPouttemp2;

            ERPs = EP_tar_post(TimeWindow{comp},chdraw_use,trial_use,s);
            [ERPaligned, Lags,ERPtrhistory,Lagtrhistory,CORs] = Woodys(ERPs,limits,threshold,Fs,Template,ERPout);

            STDs(comp,end,s) = std(Lags(:)/Fs*1000);
            LAGs(trial_use,comp,end,s) = Lags(:);


            ERPsLong = EP_tar_post(:,chdraw_use,trial_use,s);
            ERPalignedLong = [];
            for tt = 1:length(trial_use)
                L = LAGs(tt,comp,end,s);
                if L < 0
                    ERPalignedLong(:,:,tt) = cat(1,ERPsLong(abs(L)+1:end,:,tt), ERPsLong(1:abs(L),:,tt));
                elseif L>=0
                    ERPalignedLong(:,:,tt) = cat(1,ERPsLong(end-abs(L)+1:end,:,tt), ERPsLong(1:end-abs(L),:,tt));
                end

            end
            ERPalsingle{comp}(:,trial_use,ismember(chdraw,chdraw_use),end,s) = permute(ERPalignedLong,[1 3 2]);
            ERPorsingle{comp}(:,trial_use,ismember(chdraw,chdraw_use),end,s) = permute(ERPsLong,[1 3 2]);
            ERPalmean{comp}(:,ismember(chdraw,chdraw_use),end,s) = mean(ERPalignedLong,3);
            ERPormean{comp}(:,ismember(chdraw,chdraw_use),end,s) = mean(ERPsLong,3);
            end
        
    end
end
%%
save(['ERP\value_LatencyShift.mat'],'LAGs')
save(['ERP\value_LatencyShiftedERP.mat'],'ERPalsingle')
save(['ERP\value_LatencyShifted_originalERP.mat'],'ERPorsingle')

%%
load([preptype,'\value_LatencyShift.mat'])
load([preptype,'\value_LatencyShiftedERP.mat'])
load([preptype,'\value_LatencyShifted_originalERP.mat'])

STDs = squeeze(std(LAGs,'omitnan'));
for comp = 1:Ncomp
    ERPalmean{comp} = squeeze(mean(ERPalsingle{comp},2,'omitnan'));

end
%%
%%

color2 = {[0 0 0];[0, 0, 255]./255;[255, 0, 0]./255;[0, 128, 0]./255;[255, 165, 0]./255;[128, 0, 128]./255;[0 0 0]};
CondName = ['Pre',CondName2(2:end),'Post'];
%-- LV
p = []; tbl = []; stats = [];
ax = [];
figure;
for comp= 1:Ncomp
    fig{comp} = subplot(1,Ncomp,comp);
    hold on;
    for sess = 1:Nsess
        ax{comp,sess} = boxchart(ones(Nsub,1)*sess,squeeze(STDs(comp,sess,:)),...
            'BoxFaceColor',color2{sess},'MarkerColor',color2{sess});

        set(gca,'FontSize',15,'xtick',1:Nsess,'xticklabel',CondName,'TickDir','out');
        xlim([0.5 7.5])
        ylabel('Latency variability')
        title(CompName{comp})
    end

    Xs = squeeze(STDs(comp,:,:))';
    for ii = 1:size(Xs,2)
       rem = isnan(Xs(:,ii));
       Xs(rem,:) = [];
    end
    [p(comp),tbl{comp},stats{comp}] = friedman(Xs,1,'off');

    post{comp} = multcompare(stats{comp},"CType","dunn-sidak",'Display','off');

end
mean(STDs,3)
for ch = 1:3; chi(ch) = tbl{ch}{2,5}; end




interval = 0.06;
for comp = 1:Ncomp

    drawsignificance(fig{comp},squeeze(STDs(comp,ConditionOrder,:))',post{comp},interval)
end
set(gcf,'Position',[200 400 950 270])

%-- LV histogram
figure;
for sess = 1:Ncon
    for comp = 1:Ncomp
        subplot(3,Ncon,(comp-1)*Ncon+sess)
        X = squeeze(LAGs(:,comp,ConditionOrder(sess),:)/Fs*1000); % in ms

        histogram(X(:),20,"Normalization","pdf",'FaceColor',color_use{ConditionOrder(sess)})
        ylim([0 0.015])

        set(gca,'fontsize',8)
        if sess == 1
            xlabel('Latency shift (ms)')
            ylabel('Probability density')
        end
        if comp == 1
            title({CondName{ConditionOrder(sess)},''})
        end

    end
end

f_ = [];
figure;
for comp = 1:Ncomp
    subplot(1,Ncomp,comp);
    hold on;
    for sess = 1:Ncon
        [f_(:,comp,sess),xi]= ksdensity(squeeze(STDs(comp,ConditionOrder(sess),:)));


        plot(xi,squeeze(f_(:,comp,sess)),'color',color_use{ConditionOrder(sess)},'linewidth',2)


    end


end

%% LV vs. amplitude
load([preptype,'\value_ERPComp_FC1FC2.mat']);
r_ = []; p_ = [];
figure;
for sess = 1:Ncon
    for comp = 1:Ncomp

        if comp == 1
            Y = squeeze(mean(P2amp(:,:,ConditionOrder(sess)),'omitnan'))';
        elseif comp == 2
            Y = squeeze(mean(P3amp(:,:,ConditionOrder(sess)),'omitnan'))';
        elseif comp == 3
            Y = squeeze(mean(-N2amp(:,:,ConditionOrder(sess)),'omitnan'))';
        end


        X= squeeze(STDs(comp,ConditionOrder(sess),:));
        [r_(comp,sess),p_(comp,sess)] = corr(X,Y,'rows','pairwise','type','Spearman');

        subplot(3,Ncon,(comp-1)*Ncon+sess)

        scatter(X,Y,[],color_use{ConditionOrder(sess)},'filled','o',...
            'MarkerEdgeColor',color_use{ConditionOrder(sess)}+0.08,'MarkerFaceAlpha',0.8);
        set(gca,'FontSize',15,'TickDir','out');
        if sess == 1
            ylabel({CompName{comp},'Amplitude '})

            xlabel('Latency variability')

        end
        if comp == 1
            title(CondName(ConditionOrder(sess)))
            %             xlim([45 85])
        elseif comp == 2
            %             xlim([45 85])

        elseif comp ==3
            %             xlim([30 85])
            %             xticks(30:20:80)

        end

    end

end
set(gcf,'Position',[680   284   870   694])


%% LV vs. Acc
load([preptype,'\value_Acc.mat'])
r_2 = []; p_2 = [];
figure;
for sess = 1:Ncon
    for comp = 1:Ncomp

        X= squeeze(STDs(comp,ConditionOrder(sess),:));
        Y = Acc_percent(:,sess);
        [r_2(comp,sess),p_2(comp,sess)] = corr(X,Y,'rows','pairwise','type','Spearman');

        subplot(3,Ncon,(comp-1)*Ncon+sess)

        scatter(X,Y,[],color_use{ConditionOrder(sess)},'filled','o',...
            'MarkerEdgeColor',color_use{ConditionOrder(sess)}+0.08,'MarkerFaceAlpha',0.8);
        set(gca,'FontSize',15,'TickDir','out','YLim',[40 105]);
        if sess == 1
            ylabel({CompName{comp},'Accuracy (%)'})

            xlabel('Latency variability')

        end
        if comp == 1
            title(CondName(ConditionOrder(sess)))

        end

    end

end
set(gcf,'Position',[680   284   870   694])


%% mediation analysis
X1 = [squeeze(STDs(3,4,:)) abs(mean(N2amp(:,:,4),1,'omitnan'))' Acc_percent(:,1) ];
save('value_Control_N2ampLV_Acc.mat','X1')

X2 = [squeeze(STDs(3,1,:)) abs(mean(N2amp(:,:,1),1,'omitnan'))' Acc_percent(:,2) ];
save('value_Listen_N2ampLV_Acc.mat','X2')

X3 = [squeeze(STDs(3,2,:)) abs(mean(N2amp(:,:,2),1,'omitnan'))' Acc_percent(:,3) ];
save('value_Speak_N2ampLV_Acc.mat','X3')
