%% Data separability
LAMBDA = 0.02;
%%
fpath = 'E:\1_연구\분석코드\2024_adaptive_BCI\Data\Pilot';
codepath = 'E:\1_연구\분석코드\2024_adaptive_BCI\AdaptiveBCI_AR\0_online_code';

SubNameListAll = {'Subtest02','Subtest03','Subtest04','Subtest05','Subtest06','Subtest07',...
    'Subtest08','Subtest09','Subtest10','Subtest11','Subtest12','Subtest13','Subtest14','Subtest15'};
remsub = find(ismember(SubNameListAll,{'Subtest06','Subtest11'}));
SubNameList = SubNameListAll(setdiff(1:length(SubNameListAll),remsub));
% SubNameList = {'Subtest02','Subtest03','Subtest04','Subtest05','Subtest07',...
%     'Subtest08','Subtest09','Subtest10','Subtest12','Subtest13','Subtest14','Subtest15'};


Nsub = length(SubNameList);

Ntr_tr = 15;
Ntr_con = 15;
Nsess = 7;
Ntr_te = Ntr_con*Nsess;

Nch = 31;
Fs = 500;
Times = -0.2:1/Fs:0.6-1/Fs;
chanloc = readlocs('D:\[1]EEGBCI\[2]Research\Code&Algorithm\BP_channelLocs.locs');

% Define the start and end colors
start_color = [0, 0, 255]; % Blue
end_color = [255, 0, 0]; % Red

gradient = zeros(Nsub, 3); % Initialize the matrix for RGB values
for i = 1:3
    gradient(:, i) = linspace(start_color(i), end_color(i), Nsub);
end
colors = gradient./255;

Nstim = 4;
Nrepeat = 10;
Scores_tar = NaN(Nrepeat,Ntr_con*Nsess,Nsub);
Scores_nar = NaN(Nrepeat,Nstim-1,Ntr_con*Nsess,Nsub);
Targets = NaN(Ntr_con*Nsess,Nsub);
Prediction = NaN(Ntr_con*Nsess,Nsub);
Updated  = NaN(Ntr_con*Nsess,Nsub);
Nupdated  = NaN(Ntr_con*Nsess,Nsub);
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

load("TrialUsed.mat")
% Trials_ = Trials(setdiff(1:length(SubNameListAll),remsub));
%%
for s = 9:Nsub
    SubName = SubNameList{s};

    %% load (Epoch)

    %% svm classification score
    load([fpath,'\Simulate_',num2str(LAMBDA),'param\',SubName])
    load([fpath,'\',SubName])
    p = load([fpath,'\Dat_',SubName,'\param.mat']);

    uptrial_id = find(Trials_{s} == Ntr_con)+1:length(Trials_{s});
    uptrial = Trials_{s}(find(Trials_{s} == Ntr_con)+1:length(Trials_{s}));

    targetall = vars.target(:);
    targetlist = targetall(uptrial);

    
    if length(UPmodel) > (Nsess-1)*Ntr_con
        DSP_init = UPmodel(1).DSP_init;
        MDL_init = p.param.trD.mdl_init;
        UPmodel(1:find(Trials_{s} == Ntr_con)) = []; % remove presession
    end

    for t = 1:length(UPmodel)

        score = UPmodel(t).classprob;
        if ~isfield(UPmodel,'target')
            target = targetlist(t);
        else
            target = UPmodel(t).target;
        end
        Scores_tar(:,uptrial(t),s) = score(:,target);
        Scores_nar(:,:,uptrial(t),s) = score(:,setdiff(1:Nstim,target));

        Targets(uptrial(t),s) = target;
        Prediction(uptrial(t),s) = UPmodel(t).prediction;
        Updated(uptrial(t),s) = UPmodel(t).upcheck;
        Nupdated(uptrial(t),s) = length(UPmodel(t).upids);
    end

end

%% d_prime
Scores_tar_ = Scores_tar(:,Ntr_con+1:end,:);
Scores_nar_ = Scores_nar(:,:,Ntr_con+1:end,:);

Correct_tar = Scores_tar_>0.5;
Corr_tar = reshape(Correct_tar,Nrepeat,Ntr_con,Nsess-1,Nsub);
Correct_nar = Scores_nar_<0.5;
Corr_nar = reshape(Correct_nar,Nrepeat*3,Ntr_con,Nsess-1,Nsub);

rec = squeeze(sum(Corr_tar)./(sum(Corr_tar)+sum(~Corr_tar)));


H_ = rec;
F_ = squeeze(sum(~Corr_nar)./(sum(Corr_nar)+sum(~Corr_nar)));

d_prime = calculateDPrime(H_,F_);

%% f1 score
sublist = 9:Nsub;

TP = sum(Corr_tar);
FP = sum(~Corr_nar);
FN = sum(~Corr_tar);
TN = sum(Corr_nar);

Precision = TP./(TP+FP) + 0.000000001;
Recall = TP./(TP+FN) + 0.000000001;
f1 = squeeze(2.*Precision.*Recall./(Precision+Recall));

f1Mean = squeeze(mean(f1,'omitnan'));

figure; plot(reshape(f1(:,:,sublist),Ntr_con*(Nsess-1),[]))

%%
Names = {'Pilot01','Pilot02','Pilot03','Pilot04','Pilot06','Pilot07','Pilot08','Pilot09','Pilot11','Pilot12','Pilot13','Pilot14'};


Correctness = Targets == Prediction;
Correctness_ = Correctness(Ntr_con+1:end,sublist);
Updated_ = Updated(Ntr_con+1:end,sublist);
Updated_(isnan(Updated_)) = 0;
CU = Correctness_ & Updated_;
WU = ~Correctness_ & Updated_;
CN = Correctness_ & ~Updated_;
WN = ~Correctness_ & ~Updated_;

f1_ = reshape(f1(:,:,sublist),[],length(sublist));
figure; 
binrange = 0:0.2:0.8;
ylimit = [0 31];
for ss = 1:2
    subplot(4,2,ss);hold on;
    histogram(f1_(CU(:,ss),ss),binrange); 
 histogram(f1_(WU(:,ss),ss),binrange); 
 ylim(ylimit)
 title([Names{sublist(ss)},' Selected'])
 subplot(4,2,ss+2);hold on;
  histogram(f1_(CN(:,ss),ss),binrange); 
   histogram(f1_(WN(:,ss),ss),binrange); 
    ylim(ylimit)
 title([Names{sublist(ss)},' Not selected'])

end
for ss = 3:4
    subplot(4,2,ss+2);hold on;
    histogram(f1_(CU(:,ss),ss),binrange); 
 histogram(f1_(WU(:,ss),ss),binrange); 
  ylim(ylimit)
 title([Names{sublist(ss)},' Selected'])

 subplot(4,2,ss+4);hold on;
  histogram(f1_(CN(:,ss),ss),binrange); 
   histogram(f1_(WN(:,ss),ss),binrange); 
    ylim(ylimit)
 title([Names{sublist(ss)},' Not Selected'])

end
xlabel('F1 score')
ylabel('Count')
legend({'Correct','Incorrect'})
% 
squeeze(sum(reshape(Nupdated,Ntr_con,[],Nsub),'omitnan'))
%%
%% plot
%% box plot
%% pre vs. post

%% main 1 ~ 4


%% statistics
%-- Friedman test

