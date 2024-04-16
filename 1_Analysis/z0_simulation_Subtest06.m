%% 
clear all; close all
%%
 fpath = 'E:\1_연구\분석코드\2024_adaptive_BCI\Data\Pilot';
codepath = 'E:\1_연구\분석코드\2024_adaptive_BCI\AdaptiveBCI_AR\0_online_code';

Ntr_tr = 8;
Ntr_con = 15;
Nsess = 7;
Ntr_te = Ntr_con*Nsess;
Fs = 500;
Delay = 0.7*Fs;

SubName = 'Subtest06';

Ntr_exp = [14 13 15 15 15 15 13];
Ntr_accum = cumsum(Ntr_exp);
Missed_trials = {2,[4, 5],0,0,0,0,[9 15]};
%% Data

filename = [fpath,'\RawData\',SubName,'.vhdr'];
EEG = pop_fileio(filename);


triglat = [EEG.event.latency];
trigtyp = {EEG.event.type};

trigEnd = find(ismember(trigtyp,'S 13'));

trig_tmp = zeros(1,EEG.pnts);

for id = 1:length(triglat)
    trig_tmp(triglat(id)) = str2double(trigtyp{id}(end-1:end));
end


% trigtyp_ = trig_tmp(trig_tmp~=0);
% triglat_ = find(trig_tmp~=0);
% trigStart_ = find(trigtyp_ == 57);
% 
% for tt = 22:length(trigStart_)
%     tt
%     [trigtyp_(trigStart_(tt)-1:trigStart_(tt)+50);...
%         diff(triglat_(trigStart_(tt)-1:trigStart_(tt)+50))./500 0]
% 
%     pause;
% end

%% training data
mkdir([fpath,'\Dat_',SubName])
for tr = 1:Ntr_tr
    % block start point
    % last end trigger (13) + (5 + 1) sample points

    if tr >1 
        ind = trigEnd(tr-1);
        startpoint = triglat(ind) + 5 + 1;

    else
       tmp = find(ismember(trigtyp,'S 15'));
       startpoint = triglat(tmp(1)) - 10*Fs; 
    end

   
    % block end point
    % end trigger (13) + 5 sample points
    ind2 = trigEnd(tr);
    endpoint = triglat(ind2) + 5;

    sig_vec = EEG.data(:,startpoint:endpoint)./0.04883;
    trigger = trig_tmp(:,startpoint:endpoint);

    trigger_re = trigger;
    trigger_re(1:Delay) = [];
    trigger_re = [trigger_re zeros(1,Delay)];

    save([fpath,'\Dat_',SubName,'\',SubName,'_Training',num2str(tr)],'sig_vec','trigger','trigger_re');
end

%% param
cd(codepath);

param = Initialization('Def');
close;

param.SubName = SubName;

param.device = 3; %('Device to control (All:0/Doorlock:1/AirConditioner:2/Lamp:3/Bluetoothspeaker:4) ');

param.NumTrTrial = Ntr_tr;

param.Stims = [1:4];
param.NumStims = 4;

param.repeat = 10;
param.prep_factor = [1:3];

param.caltime = 60*param.Fs;

param.winsize = 5;
param.DSP.type = 'nf';
param.DSP.nf = 5; % use 5 spatial filters
param.DSP.theta = 0.1;
param.DSP.window = param.Baseline + 1:param.winsize:param.Baseline+param.Epocline;


sig = []; trig = [];
for tr = 1:Ntr_tr % data reload
    load([fpath,'\Dat_',SubName,'\',SubName,'_Training',num2str(tr)]);
    sig = cat(2,sig,sig_vec);
    trig = cat(2,trig,trigger_re);
end

%-- restore params
[C, param] = P300_processing_adaptive(sig(1:param.NumChIni,:),trig, param); % process data

param.trD.mode  = 'testing'; % change mode to test


save([fpath,'\Dat_',SubName,'\',SubName,'_Training'],'sig','trig'); % save data
save([fpath,'\Dat_',SubName,'\param'],'param'); % save parameter

%% prediction from log

logdata = importdata([fpath,'\Log_',SubName,'.txt']);
teststart = find(isnan(logdata.data)) ;
predicted = zeros(Ntr_con,Nsess);

%-- pre session
log_test = logdata.data(teststart(1)+1:teststart(1)+Ntr_exp(1)*2);    
predicted(1:length(log_test)/2,1) = log_test(2:2:end);

%-- main session
 
for ss = 1:5
    log_test2 = logdata.data(teststart+Ntr_accum(ss)*2+1:teststart+Ntr_accum(ss+1)*2); 
    predicted(1:length(log_test2)/2,1+ss) = log_test2(2:2:end);
end

%-- post session
log_test3 = logdata.data(teststart+Ntr_accum(6)*2+1:teststart+Ntr_accum(7)*2);
predicted(1:length(log_test3)/2,end) = log_test3(2:2:end);

%% pre-session (1 ~ 15)
param.trD.threshold = 0.5;
param.trD.adaptmode = 'margin';
param.DSP.lambda = 0.02;
param.trD.ADmode = 'fixed';
param.prediction = [];
param.Numtrial = 0;

targetlist = [];
for tr = 1:Ntr_exp(1)
    param.Numtrial = param.Numtrial + 1;   
    % block start point
    % last end trigger (13) + (5 + 1) sample points
    ind = trigEnd(Ntr_tr+tr-1);
    startpoint = triglat(ind) + 5 + 1;

    % block end point
    % end trigger (13) + 5 sample points
    ind2 = trigEnd(Ntr_tr+tr);
    endpoint = triglat(ind2) + 5;

    sig_vec = EEG.data(:,startpoint:endpoint)./0.04883;
    trigger = trig_tmp(:,startpoint:endpoint);

    trigger_re = trigger;
    trigger_re(1:Delay) = [];
    trigger_re = [trigger_re zeros(1,Delay)];

    
   [C, param] = P300_processing_adaptive(sig_vec(1:param.NumChIni,:),trigger_re,param); % process data
    param.prediction(param.Numtrial) = C;

    save([fpath,'\Dat_',SubName,'\',SubName,'_Testing',num2str(tr)],'sig_vec','trigger','trigger_re');

    trigtype = trigger_re(trigger_re~=0);
    tartemp = trigtype(find(trigtype == 57)-1);
    targetlist(param.Numtrial,1) = tartemp(end);
end
targetlist(2,1) = NaN;
save([fpath,'\Dat_',SubName,'\param'],'param'); % save parameter


%% main session

load([fpath,'\Dat_',SubName,'\param'],'param'); % load parameter

param.trD.ADmode = 'adaptive';
param.Numtrial = Ntr_exp(1);
for sess = 2:6
for tr_s = 1:Ntr_exp(sess)
    param.Numtrial = param.Numtrial + 1;   
    % block start point
    % last end trigger (13) + (5 + 1) sample points
    ind = trigEnd(Ntr_tr+param.Numtrial-1);
    startpoint = triglat(ind) + 5 + 1;

    % block end point
    % end trigger (13) + 5 sample points
    ind2 = trigEnd(Ntr_tr+param.Numtrial);
    endpoint = triglat(ind2) + 5;

    sig_vec = EEG.data(:,startpoint:endpoint)./0.04883;
    trigger = trig_tmp(:,startpoint:endpoint);

    trigger_re = trigger;
    trigger_re(1:Delay) = [];
    trigger_re = [trigger_re zeros(1,Delay)];

   [C, param] = P300_processing_adaptive(sig_vec(1:param.NumChIni,:),trigger_re,param); % process data
    param.prediction(param.Numtrial) = C;

    save([fpath,'\Dat_',SubName,'\',SubName,'_Testing',num2str(param.Numtrial)],'sig_vec','trigger','trigger_re');

    trigtype = trigger_re(trigger_re~=0);
    tartemp = trigtype(find(trigtype == 57)-1);
    if ~isempty(tartemp)
        targetlist(tr_s,sess) = tartemp(end);
    else
        targetlist(tr_s,sess) = NaN;
    end
end
end
save([fpath,'\Dat_',SubName,'\param'],'param'); % save parameter


%% post session 
load([fpath,'\Dat_',SubName,'\param'],'param'); % load parameter
param.Numtrial = Ntr_accum(6);
for tr = 1:Ntr_exp(end)
    param.Numtrial = param.Numtrial + 1;   
    % block start point
    % last end trigger (13) + (5 + 1) sample points
    ind = trigEnd(Ntr_tr+param.Numtrial-1);
    startpoint = triglat(ind) + 5 + 1;

    % block end point
    % end trigger (13) + 5 sample points
    ind2 = trigEnd(Ntr_tr+param.Numtrial);   
    endpoint = triglat(ind2) + 5;

    sig_vec = EEG.data(:,startpoint:endpoint)./0.04883;
    trigger = trig_tmp(:,startpoint:endpoint);

    trigger_re = trigger;
    trigger_re(1:Delay) = [];
    trigger_re = [trigger_re zeros(1,Delay)];

   [C, param] = P300_processing_adaptive(sig_vec(1:param.NumChIni,:),trigger_re,param); % process data
    param.prediction(param.Numtrial) = C;

    save([fpath,'\Dat_',SubName,'\',SubName,'_Testing',num2str(param.Numtrial)],'sig_vec','trigger','trigger_re');


     trigtype = trigger_re(trigger_re~=0);
     if tr == 8
         targetlist(tr,7) = trigtype(1);
     elseif tr == 13
         targetlist(tr,7) = NaN;
     else
    targetlist(tr,7) = trigtype(find(trigtype == 57)-1);
     end
end

save([fpath,'\Dat_',SubName,'\param'],'param'); % save parameter
%% accruacy

targetlist(targetlist == 0) = NaN;
Acc_online = sum(targetlist == predicted)./Ntr_exp;

restored_predict = NaN(Ntr_con,Nsess);
for tt = 1:length(Ntr_exp)
    if tt == 1
        ind = 1:Ntr_exp(tt);
    else
        ind = Ntr_accum(tt-1)+1:Ntr_accum(tt-1)+Ntr_exp(tt);
    end

restored_predict(1:Ntr_exp(tt),tt) = param.prediction(ind);
end

Acc_restored = sum(targetlist == restored_predict)./Ntr_exp;

vars.target = targetlist;
vars.SubName = SubName;
save([fpath,'\',SubName])
%% post session (fixed)
% param.trD.ADmode = 'fixed';
% 
% for tr = Ntr_con*(Nsess-1)+1:Ntr_con*Nsess
%     param.Numtrial = param.Numtrial + 1;
% 
%     datapath = [fpath,'\Dat_',SubName,'\',SubName,'_Testing',num2str(tr)];
%     load(datapath)
% 
%     [C, param] = P300_processing_adaptive(sig_vec(1:param.NumChIni,:),trigger_re,param);
%     
%     param.prediction(param.Numtrial) = C;
% end