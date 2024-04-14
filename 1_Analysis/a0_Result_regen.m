%% Data arrange


clear all; close all
%% parameters
fpath = 'E:\1_연구\분석코드\2024_adaptive_BCI\Data\Pilot';
codepath = 'E:\1_연구\분석코드\2024_adaptive_BCI\AdaptiveBCI_AR\0_online_code';

Ntr_tr = 8;
Ntr_con = 15;
Nsess = 7;
Ntr_te = Ntr_con*Nsess;

B_fir_lf_2 = getfirfiltcoeff(15,'low',500,0);
LPF15 = 1;

% SubNameList = {'Subtest02','Subtest03','Subtest04','Subtest05'};
SubNameList = {'Subtest06','Subtest07'};

Nsub = length(SubNameList);
%%
for s = 2:Nsub
    SubName = SubNameList{s};

    p = load([fpath,'\Dat_',SubName,'\param.mat']);


    %% % Training

    dataname = [SubName,'_Training.mat'];

    cd(codepath)

    load([fpath,'\Dat_',SubName,'\',dataname]);

    sig = double(sig*0.04883);
    [sig, param]            = PreProcess(sig,p.param);
    % 15Hz LPF
    if LPF15 == 1
        sig = filtfilt(B_fir_lf_2,1,sig')';
    end

    param.trD.mode = 'training';
    EP = EpochData(sig,trig,param);
    EP_ = Epoch_condition(EP,param);
    EP_.badch = param.badch;


    trigtype = trig(trig~=0);
    triglat = find(trig~=0);

    blockst = find(trigtype == 12);
    blockend = find(trigtype == 13);

    targets = trigtype(blockst - 1);

    Block = [];
    for b= 1:length(blockst)
        blockind = triglat(blockst(b)):triglat(blockend(b));
        sig_b = sig(:,blockind);
        trig_b = trig(:,blockind);
        Block.sig{b} = sig_b;
        Block.trig{b} = trig_b;

    end
    Block.target = targets;
    Block.badch = param.badch;

    Epoch = EP_;

    if LPF15 == 1
        save([fpath,'\Train\Block_15HZLPF\',SubName],'Block');
        save([fpath,'\Train\Epoch_15HZLPF\',SubName],'Epoch');

    else
        save([fpath,'\Train\Block\',SubName],'Block');
        save([fpath,'\Train\Epoch\',SubName],'Epoch');
    end
    %% Test
    load([fpath,'\',SubName,'.mat'])
    targets = vars.target(:);
    Block = []; Epoch = [];
    for tr = 1:Ntr_te
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
        EP.target = targets(tr);
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
    if LPF15 == 1
        save([fpath,'\Test\Block_15HZLPF\',SubName],'Block');
        save([fpath,'\Test\Epoch_15HZLPF\',SubName],'Epoch');

    else
        save([fpath,'\Test\Block\',SubName],'Block')
        save([fpath,'\Test\Epoch\',SubName],'Epoch');
    end
end

%% Subtest06 data 

load([fpath,'\Test\Block\Subtest06'])




