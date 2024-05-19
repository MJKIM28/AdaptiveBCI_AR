%% Data arrange


clear all; close all
%% parameters
analpath = cd;
fpath = 'E:\1_연구\분석코드\2024_adaptive_BCI\Data\Pilot';
codepath = 'E:\1_연구\분석코드\2024_adaptive_BCI\AdaptiveBCI_AR\0_online_code';

Ntr_tr = 15;
Ntr_con = 15;
Nsess = 7;

B_fir_lf_2 = getfirfiltcoeff(15,'low',500,0);
LPF15 = 0;

SubNameList = {'Subtest02','Subtest03','Subtest04','Subtest05',...
'Subtest07','Subtest08','Subtest09','Subtest10','Subtest12','Subtest13','Subtest14','Subtest15'}; % Ntr_con = 15

SubNameList2 = {'Subtest17','Subtest18','Subtest19','Subtest20','Subtest21','Subtest22'}; % Ntr_con = 14

SubListUse = SubNameList;
Nsub = length(SubListUse);


load('TrialUsed.mat')
trial_use = Trials;

cd(codepath)
%%
for s = 1:Nsub
if s <5
    Nsess = 6;
else
    Nsess = 7;
end

    SubName = SubListUse{s};

    p = load([fpath,'\Dat_',SubName,'\param.mat']);


    %% % Training

%     dataname = [SubName,'_Training.mat'];
% 

% 
%     load([fpath,'\Dat_',SubName,'\',dataname]);
% 
%     sig = double(sig*0.04883);
%     [sig, param]            = PreProcess(sig,p.param);
%     % 15Hz LPF
%     if LPF15 == 1
%         sig = filtfilt(B_fir_lf_2,1,sig')';
%     end
% 
%     param.trD.mode = 'training';
%     EP = EpochData(sig,trig,param);
%     EP_ = Epoch_condition(EP,param);
%     EP_.badch = param.badch;
% 
% 
%     trigtype = trig(trig~=0);
%     triglat = find(trig~=0);
% 
%     blockst = find(trigtype == 12);
%     blockend = find(trigtype == 13);
% 
%     targets = trigtype(blockst - 1);
% 
%     Block = [];
%     for b= 1:length(blockst)
%         blockind = triglat(blockst(b)):triglat(blockend(b));
%         sig_b = sig(:,blockind);
%         trig_b = trig(:,blockind);
%         Block.sig{b} = sig_b;
%         Block.trig{b} = trig_b;
% 
%     end
%     Block.target = targets;
%     Block.badch = param.badch;
% 
%     Epoch = EP_;
% 
%     if LPF15 == 1
%         save([fpath,'\Train\Block_15HZLPF\',SubName],'Block');
%         save([fpath,'\Train\Epoch_15HZLPF\',SubName],'Epoch');
% 
%     else
%         save([fpath,'\Train\Block\',SubName],'Block');
%         save([fpath,'\Train\Epoch\',SubName],'Epoch');
%     end
    %% Test
    load([fpath,'\',SubName,'.mat'])
    targets = vars.target(:);
    Block = []; Epoch = [];
    for tr = 1:length(trial_use{s})
        dataname =  [SubName,'_Testing',num2str(tr),'.mat'];
        load([fpath,'\Dat_',SubName,'\',dataname]);

        sig = double(sig_vec*0.04883);
        [sig, param]            = PreProcess(sig,p.param);
        % 15Hz LPF
        if LPF15 == 1
            sig = filtfilt(B_fir_lf_2,1,sig')';
        end

        EP = EpochData(sig,trigger_re,param);

        param.trD.mode = 'training';
        EP.target = targets(trial_use{s}(tr));
        EP_ = Epoch_condition(EP,param);


        trigtype = trigger_re(trigger_re~=0);
        triglat = find(trigger_re~=0);

        blockst = find(trigtype == 12);
        blockend = find(trigtype == 13);


        blockind = triglat(blockst):triglat(blockend);
        sig_b = sig(:,blockind);
        trig_b = trigger_re(:,blockind);

        Block.sig{tr} = sig_b;
        Block.trig{tr} = trig_b;

        Epoch{tr} = EP_;

    end

    Block.target = targets;
    Block.badch = param.badch;
    Epoch{1}.badch = param.badch;
    if LPF15 == 1
        save([fpath,'\Test\Block_15HZLPF\',SubName],'Block');
        save([fpath,'\Test\Epoch_15HZLPF\',SubName],'Epoch');

    else
        save([fpath,'\Test\Block\',SubName],'Block')
        save([fpath,'\Test\Epoch\',SubName],'Epoch');
    end
end
cd(analpath);


%%  save missed trial
% Subtest02 ~ Subtest16
% Subtest06, Subtest11, Subtest16 제외
Trials = cell(12,1); 
for s = 1:12
    if s < 5
        Nsess = 6;
    else
        Nsess = 7;
    end
    Trials{s} = 1:Ntr_con*Nsess;
end
% for s = 1:4
%     Trials{s} = [1:Ntr_con*5, Ntr_con*(Nsess-1)+1:Ntr_con*Nsess]; % Subtest02 ~ Subtest05: 4 main session 
% end
Trials{6} = [1:4,7:Ntr_con*Nsess];% Subtest08?? pre 5,6th block missed 
Trials{11} = [1:15,17:Ntr_con*Nsess]; % Subtest14 main 1 1st block missed

save('TrialUsed','Trials')

%%  save missed trial ver2 (14 trials)
% % Subtest17 ~ Subtest22
% Ntr_con2 = 14;
% Trials2 = cell(6,1); 
% for s = 1:6
%     Trials2{s} = 1:Ntr_con2*Nsess;
% end
% Trials2{1} = [1:14, 29:Ntr_con2*7]; % Subtest17 main 2 ~ (main1 x)
% Trials2{3} = [1:Ntr_con2*6]; % Subtest19 post x
% Trials2{4} = [1:Ntr_con2*6]; % Subtest20 post x
% save('TrialUsed2','Trials2')


