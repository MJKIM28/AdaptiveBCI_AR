%% alpha modulation

%% parameters
fpath = 'E:\1_연구\분석코드\2024_adaptive_BCI\Data\Pilot';
codepath = 'E:\1_연구\분석코드\2024_adaptive_BCI\AdaptiveBCI_AR\0_online_code';


SubNameList = {'Subtest02','Subtest03','Subtest04','Subtest05','Subtest07',...
    'Subtest08','Subtest09','Subtest10','Subtest12','Subtest13','Subtest14','Subtest15'};
SubNameList2 = {'Subtest17','Subtest18','Subtest19','Subtest20','Subtest21','Subtest22'};

Nsub1 = length(SubNameList);
Nsub2 = length(SubNameList2);

Ntr_con = 14;
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


Bands = [[1 4];[4 8]; [8 12];[12 30]];
BandName = {'Delta','Theta','Alpha','Beta'};
Nband = length(Bands);

SessName = {'Pre','Main 1','Main 2','Main 3','Main 4','Main 5','Post (adaptive)'};

%% alpha power computation
%-- all channel, all condition
n12 = load(['value_Power_post_',num2str(Nsub1),'subs.mat']);
n6 = load(['value_Power_post_',num2str(Nsub2),'subs.mat']);
f = n6.f;
Ptemp1 = n12.Ppost(:,:,:,:,9:end); % subtest12 ~ subtest15
Ptemp2 = cat(3,NaN(Nf,Nch,1,Nsess,Nsub2),n6.Ppost);

Ppost = cat(5,Ptemp1,Ptemp2);
Nsub = size(Ppost,5);
P_alpha1_d = squeeze(mean(trapz(df,Ppost(f>Bands(3,1) & f<Bands(3,2),[5 6 7],:,:,:)),2,'omitnan'));
P_alpha2_d = squeeze(mean(trapz(df,Ppost(f>Bands(3,1) & f<Bands(3,2),[25 26 27],:,:,:)),2,'omitnan'));
P_alpha3_d = squeeze(mean(trapz(df,Ppost(f>Bands(3,1) & f<Bands(3,2),[29 30 31],:,:,:)),2,'omitnan'));

PalphaBCI{1} = squeeze(mean(P_alpha1_d(:,1,:),'omitnan'))';
PalphaBCI{2} = squeeze(mean(P_alpha2_d(:,1,:),'omitnan'))';
PalphaBCI{3} = squeeze(mean(P_alpha3_d(:,1,:),'omitnan'))';


%-- alpha power modulation (compared to pre session)
Pratio = [];
Pratio{1} = squeeze(10*log10(P_alpha1_d./reshape(PalphaBCI{1},1,1,Nsub)));
Pratio{2} = squeeze(10*log10(P_alpha2_d./reshape(PalphaBCI{2},1,1,Nsub)));
Pratio{3} = squeeze(10*log10(P_alpha3_d./reshape(PalphaBCI{3},1,1,Nsub)));

%%
%% plot
%% box plot 


APM = squeeze(mean(Pratio{3},1,'omitnan'))';
sub_sp = 1:10;
% sub_wat = 11:Nsub;
figure;hold on;
        a1 = [];
        for con = 1:Nsess
            a1(con) = boxchart((con)*ones(length(sub_sp),1),APM(sub_sp,con),...
                'BoxWidth',0.2,'MarkerStyle','.');
        end

    set(gca,'xtick',1:Nsess,'xticklabel',SessName,'fontsize',15,'TickDir','out')
    ylim([min(min(APM))*1.2 max(max(APM))*1.2])

    ylabel('dB')

    figure;hold on;
        a1 = [];
        for con = 2:Nsess-1
            a1(con) = boxchart((con)*ones(length(sub_wat),1),APM(sub_wat,con),...
                'BoxWidth',0.2,'MarkerStyle','.');
        end

    set(gca,'xtick',1:Nsess,'xticklabel',SessName,'fontsize',15,'TickDir','out')
    ylim([min(min(APM))*1.2 max(max(APM))*1.2])

    ylabel('dB')


%% pre vs. post

%% main 1 ~ 4

%% control vs DE

%% statistics
%-- Friedman test

FM = [];
for ch = 1:3
 xx= Pratio{ch};
 isnan(xx)
    [p,tbl,stats] = friedman(,1,"off");
    FM{ch} = tbl;
    PH{ch} = multcompare(stats,"Display","off",'CType','dunn-sidak')

end
