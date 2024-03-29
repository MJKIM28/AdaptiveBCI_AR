%% Performance re-generation



%% parameters
fpath = 'E:\1_연구\분석코드\2024_adaptive_BCI\Data\Pilot';
codepath = 'E:\1_연구\분석코드\2024_adaptive_BCI\AdaptiveBCI_AR\0_online_code';

Ntr_tr = 8;
Ntr_con = 15;
Nsess = 6;
Ntr_te = Ntr_con*Nsess;

SubName = 'Subtest02';



%% Training 
p = load([fpath,'\Dat_',SubName,'\param.mat']);

dataname = [SubName,'_Training.mat'];

cd(codepath)

load([fpath,'\Dat_',SubName,'\',dataname]);

sig = double(sig*0.04883);
[sig, param]            = PreProcess(sig,p.param);

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


save([fpath,'\Train\Block\',SubName],'Block');
save([fpath,'\Train\Epoch\',SubName],'Epoch');

%% Test
load([fpath,'\',SubName,'.mat'])
targets = vars.target(:);
Block = []; Epoch = [];
for tr = 1:Ntr_te
    dataname =  [SubName,'_Testing',num2str(tr),'.mat'];
    load([fpath,'\Dat_',SubName,'\',dataname]);

    sig = double(sig_vec*0.04883);
    [sig, param]            = PreProcess(sig,p.param);


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

save([fpath,'\Test\Block\',SubName],'Block')
save([fpath,'\Test\Epoch\',SubName],'Epoch');
    
%% Accuracy 
%% fixed 
%-- train
load([fpath,'\',SubName,'.mat'])
targets = vars.target(:);

e = load([fpath,'\Train\Epoch\',SubName]);

[Feature, label, param] = FeatureExtraction(e.Epoch, param);
[C,param_in]               = Classification(Feature,label, param);


%-- test


et = load([fpath,'\Test\Epoch\',SubName]);
pred = [];
for tr = 1:Ntr_te


    [Feature, label, param_in] = FeatureExtraction(et.Epoch{tr}, param_in);
    [C,param_in]               = Classification(Feature,label, param_in);
    pred(tr) = C;

end

targets_con = reshape(targets,[Ntr_con,Nsess]);
pred_con = reshape(pred,[Ntr_con,Nsess]);

Acc = [];
for con = 1:Nsess
Acc(con) = mean(targets_con(:,con) == pred_con(:,con));

end

%% adpative 
pred_ad = [];
for tr = 1:Ntr_te


    paramNew = p.param;
    if tr > 1
        if isfield(p.param.update{tr-1},'DSP')
            paramNew.DSP = p.param.update{tr-1}.DSP;
            paramNew.trD.mdl = p.param.update{tr-1}.mdl.mdl;
        else
            paramNew.DSP = p.param.update{tr-2}.DSP;
            paramNew.trD.mdl = p.param.update{tr-2}.mdl.mdl;
        end
    end
    [Feature, label, paramNew] = FeatureExtraction(et.Epoch{tr}, paramNew);
    C               = Classification(Feature,label, paramNew);
    pred_ad(tr) = C;

end

pred_ad_con = reshape(pred_ad,[Ntr_con,Nsess]);

Acc_ad = [];
for con = 1:Nsess
    Acc_ad(con) = mean(targets_con(:,con) == pred_ad_con(:,con));
end


%% Simulation
