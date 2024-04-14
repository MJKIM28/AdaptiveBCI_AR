%% 
clear all; close all
%%
 fpath = 'E:\1_연구\분석코드\2024_adaptive_BCI\Data\Pilot';
codepath = 'E:\1_연구\분석코드\2024_adaptive_BCI\AdaptiveBCI_AR\0_online_code';

Ntr_tr = 8;
Ntr_con = 15;
Nsess = 6;
Ntr_te = Ntr_con*Nsess;
Fs = 500;
Delay = 0.7*Fs;

SubName = 'Subtest03';


%% Data
% 
% filename = [fpath,'\RawData\',SubName,'.vhdr'];
% EEG = pop_fileio(filename);
% 
% 
% triglat = [EEG.event.latency];
% trigtyp = {EEG.event.type};
% 
% trigEnd = find(ismember(trigtyp,'S 13'));
% 
% trig_tmp = zeros(1,EEG.pnts);
% 
% for id = 1:length(triglat)
%     trig_tmp(triglat(id)) = str2double(trigtyp{id}(end-1:end));
% end
%% pre-session (1 ~ 15), main 1 (16 ~30)
% for tr = 1:30
%     % block start point
%     % last end trigger (13) + (5 + 1) sample points
%     ind = trigEnd(Ntr_tr+tr-1);
%     startpoint = triglat(ind) + 5 + 1;
% 
%     % block end point
%     % end trigger (13) + 5 sample points
%     ind2 = trigEnd(Ntr_tr+tr);
%     endpoint = triglat(ind2) + 5;
% 
%     sig_vec = EEG.data(:,startpoint:endpoint)./0.04883;
%     trigger = trig_tmp(:,startpoint:endpoint);
% 
%     trigger_re = trigger;
%     trigger_re(1:Delay) = [];
%     trigger_re = [trigger_re zeros(1,Delay)];
% 
%     save(['Dat_',SubName,'\',SubName,'_Testing',num2str(tr)],'sig_vec','trigger','trigger_re');
% end

%%
%% Simulation
 p = load([fpath,'\Dat_',SubName,'\param.mat']);
 close;
param = p.param;
param.trD.mdl = p.param.trD.mdl_init;
param.trD.ADmode = 'fixed';
param.trD.feature = p.param.trD.feature(1:40*8,:);
param.trD.label = p.param.trD.label(1:40*8);
param.DSP = p.param.DSP.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old.old;
%% pre session
cd(codepath);

param.Numtrial = 0;
for tr = 1:Ntr_con
    param.Numtrial = param.Numtrial + 1;

    datapath = [fpath,'\Dat_',SubName,'\',SubName,'_Testing',num2str(tr)];
    load(datapath)

    [C, param] = P300_processing_adaptive(sig_vec(1:param.NumChIni,:),trigger_re,param);
    
    param.prediction(param.Numtrial) = C;
end

%% main 1 ~ main 4
param.trD.ADmode = 'adaptive';

for tr = Ntr_con+1:Ntr_con*(Nsess-1)
    param.Numtrial = param.Numtrial + 1;

    datapath = [fpath,'\Dat_',SubName,'\',SubName,'_Testing',num2str(tr)];
    load(datapath)

    [C, param] = P300_processing_adaptive(sig_vec(1:param.NumChIni,:),trigger_re,param);
    
    param.prediction(param.Numtrial) = C;
end


%% post session (adaptive)
param.trD.ADmode = 'adaptive';

for tr = Ntr_con*(Nsess-1)+1:Ntr_con*Nsess
    param.Numtrial = param.Numtrial + 1;

    datapath = [fpath,'\Dat_',SubName,'\',SubName,'_Testing',num2str(tr)];
    load(datapath)

    [C, param] = P300_processing_adaptive(sig_vec(1:param.NumChIni,:),trigger_re,param);
    
    param.prediction(param.Numtrial) = C;
end

save([fpath,'\Dat_',SubName,'\param.mat'],'param');
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