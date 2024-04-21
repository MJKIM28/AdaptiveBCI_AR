subname = 'Subtest11';
Folder = ['Dat_',num2str(subname)];
% load([Folder,'/param.mat'])
param = Initialization('Def');
param.SubName = subname;
param.winsize = 5;
param.DSP.type = 'nf';
param.DSP.nf = 5; % use 5 spatial filters
param.DSP.theta = 0.1;
param.DSP.window = param.Baseline + 1:param.winsize:param.Baseline+param.Epocline;

%%
% param.DSP = rmfield(param.DSP,'W');
% param.trD = rmfield(param.trD,'mdl');
% param.trD.mode = 'training';


sig = []; trig = [];
for n = 1:15 % data reload
    load([Folder,'/',param.SubName,'_Training',num2str(n)]);
    sig = cat(2,sig,sig_vec);
    trig = cat(2,trig,trigger_re);
end
clear sig_vec trigger
[C, param] = P300_processing_adaptive_simul(sig(1:param.NumChIni,:),trig, param); % process data

%%

load([subname])

param.trD.ADmode = 'fixed';
prediction0 = [];
tt = 1;
for n = 1:15
     load([Folder,'/',param.SubName,'_Testing',num2str(n)]);
   [C, param] = P300_processing_adaptive_simul(sig_vec(1:param.NumChIni,:),trigger_re,param); % process data

   prediction0(tt) = C;
    tt = tt+1;
end

prediction_00 = reshape(prediction0,15,[]);
Acc_pre = mean(vars.target(:,1) == prediction_00)

%%
load([subname])

param.trD.ADmode = 'fixed';
prediction0 = [];
tt = 1;
for n = 16:105
    param.Numtrial = param.Numtrial +1;
     load([Folder,'/',param.SubName,'_Testing',num2str(n)]);
   [C, param] = P300_processing_adaptive_simul(sig_vec(1:param.NumChIni,:),trigger_re,param); % process data

   prediction0(tt) = C;
    tt = tt+1;
end

prediction_0 = reshape(prediction0,15,[]);
Acc_fix = mean(vars.target(:,2:end) == prediction_0)


%%
param.trD.ADmode = 'adaptive';
param.trD.threshold = 0.5;
param.trD.adaptmode = 'margin';
param.DSP.lambda = 0.02;

prediction = [];
tt = 1;
param.Numtrial = 15;
for n = 16:105
    param.Numtrial = param.Numtrial +1;
     load([Folder,'/',param.SubName,'_Testing',num2str(n)]);
   [C, param] = P300_processing_adaptive_simul(sig_vec(1:param.NumChIni,:),trigger_re,param); % process data

   prediction(tt) = C;
    tt = tt+1;
end

prediction_ = reshape(prediction,15,[]);
Acc_ad = mean(vars.target(:,2:end) == prediction_)

