%% Data separability
LAMBDA = 0.02;
%%
fpath = 'E:\1_연구\분석코드\2024_adaptive_BCI\Data\Pilot';
codepath = 'E:\1_연구\분석코드\2024_adaptive_BCI\AdaptiveBCI_AR\0_online_code';

SubNameList = {'Subtest02','Subtest03','Subtest04','Subtest05','Subtest06','Subtest07'};

Nsub = length(SubNameList);

Ntr_tr = 8;
Ntr_con = 15;
Nsess = 6;
Ntr_te = Ntr_con*Nsess;

Nch = 31;
Fs = 500;
Times = -0.2:1/Fs:0.6-1/Fs;
chanloc = readlocs('D:\[1]EEGBCI\[2]Research\Code&Algorithm\BP_channelLocs.locs');



Nstim = 4;
Nrepeat = 10;
Scores_tar = NaN(Nrepeat,Ntr_con*6,Nsub);
Scores_nar = NaN(Nrepeat,Nstim-1,Ntr_con*6,Nsub);
Targets = NaN(Ntr_con*6,Nsub);


%%
param.DSP.type = 'nf';
param.trD.mode = 'testing';
param.repeat = 10;
param.Stims = 1:4;
param.NumStims = length(param.Stims);
param.Fs = Fs;
param.Baseline = 0.2*param.Fs;
param.Epocline = 0.6*param.Fs;
param.Totalepoc = param.Baseline + param.Epocline;
%%
for s = setdiff(1:Nsub,5)
    SubName = SubNameList{s};

    %% load (Epoch)

    %% svm classification score
    load([fpath,'\Simulate_',num2str(LAMBDA),'param\',SubName])

    for t = 1:length(UPmodel)
    score = UPmodel(t).classprob;
    target = UPmodel(t).target;

    Scores_tar(:,t,s) = score(:,target);
    Scores_nar(:,:,t,s) = score(:,setdiff(1:Nstim,target));

    Targets(t,s) = target;
    end

end

%% d_prime
Correct_tar = Scores_tar>0.5;
Corr_tar = reshape(Correct_tar,Nrepeat,Ntr_con,6,Nsub);
Correct_nar = Scores_nar<0.5;
Corr_nar = reshape(Correct_nar,Nrepeat*3,Ntr_con,6,Nsub);

rec = squeeze(sum(Corr_tar)./(sum(Corr_tar)+sum(~Corr_tar)));


H_ = rec;
F_ = squeeze(sum(~Corr_nar)./(sum(Corr_nar)+sum(~Corr_nar)));

d_prime = calculateDPrime(H_,F_);

%% f1 score
TP = sum(Corr_tar);
FP = sum(~Corr_nar);
FN = sum(~Corr_tar);
TN = sum(Corr_nar);

Precision = TP./(TP+FP) + 0.000000001;
Recall = TP./(TP+FN) + 0.000000001;
f1 = squeeze(2.*Precision.*Recall./(Precision+Recall));

f1Mean = squeeze(mean(f1,'omitnan'));
%%
%% plot
%% box plot
%% pre vs. post

%% main 1 ~ 4


%% statistics
%-- Friedman test

