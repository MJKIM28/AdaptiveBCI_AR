clear all; close all;
%%  power
fpath = 'E:\1_연구\분석코드\2024_adaptive_BCI\Data\Pilot';
codepath = 'E:\1_연구\분석코드\2024_adaptive_BCI\AdaptiveBCI_AR\0_online_code';


SubNameList = {'Subtest02','Subtest03','Subtest04','Subtest05','Subtest07','Subtest08','Subtest09','Subtest10','Subtest12','Subtest13','Subtest14','Subtest15'};
SubNameList2 = {'Subtest12','Subtest13','Subtest14','Subtest15','Subtest17',...
    'Subtest18','Subtest19','Subtest20','Subtest21','Subtest22'};
sublistuse = SubNameList;

Nsub = length(sublistuse);


Ntr_con = 15;
Nsess = 7;
Ntr_te = Ntr_con*Nsess;


ChNames = {'Fp1'    'Fpz'    'Fp2'    'F7'    'F3'    'Fz'    'F4'    'F8' ...
    'FT9'    'FC5'    'FC1'    'FC2'    'FC6'    'FT10'    'T7'    'C3'...
    'Cz'    'C4'    'T8'    'CP5'    'CP1'    'CP2'    'CP6'    'P7'...
    'P3'    'Pz'    'P4'    'P8'    'O1'    'Oz'    'O2'};
Nch = length(ChNames);

chanloc = readlocs('D:\[1]EEGBCI\[2]Research\Code&Algorithm\BP_channelLocs.locs');


Fs = 500;
WinLen = 1;
L = WinLen*Fs;
win = hanning(L);
Overlap = 0.5;
noverlap = L*Overlap;
nfft = 2^(nextpow2(L));
dt =  1/Fs * (WinLen*Fs - noverlap);

Nf = nfft/2+1;
df = Fs/nfft;

Nt = length(WinLen/2:dt:10.4-WinLen/2);
Nt_b = length(WinLen/2:dt:5-WinLen/2);


load('TrialUsed.mat')
load('TrialUsed2.mat')
% TrialUse = [Trials(9:end); Trials2]; %subtest12 ~ subtest22 
% TrialUse = Trials2; %subtest17 ~ subtest22;
TrialUse = Trials;
%% power spectrum

Ppost = NaN(Nf,Nch,Ntr_con*Nsess,Nsub); % dual task
%%

for s = 1:Nsub%setdiff(1:Nsub,[5,10])
    fprintf('%s\n',sublistuse{s})
    dat = load([fpath,'\Test\Block\',sublistuse{s},'.mat']);
    chlist = setdiff(1:Nch,dat.Block.badch);


   Ptmp = [];
    for t = 1:length(dat.Block.sig)
        sig =  dat.Block.sig{t};
        trig = dat.Block.trig{t};

        trigtype = trig(trig~=0);
        triglat = find(trig~=0);


        firststim_vis = triglat(find(trigtype == 12)+1);

        sig_post = sig(:,firststim_vis:firststim_vis+5.4*Fs-1);

        for ch = 1:size(sig_post,1)

            [ppost,f] = pwelch(sig_post(ch,:)',win,noverlap,nfft,Fs);
            Ptmp(:,ch,t) = ppost;

        end

        clear sig_post

    end

    trial_use = TrialUse{s};

    Ppost(:,chlist,trial_use,s) = Ptmp;


end
Ppost = reshape(Ppost,Nf,Nch,Ntr_con,Nsess,Nsub);

save(['value_Power_post_',num2str(Nsub),'subs.mat'],'Ppost','f')

