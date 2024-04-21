function param                  = Initialization(chtype)
% chtype [Def/AR/ARm]: type of channel config

%% channel config
param.Ch32 = {'FP1'    'FPZ'    'FP2'    'F7'    'F3'    'FZ'    'F4'    'F8' ...
            'FT9'    'FC5'    'FC1'    'FC2'    'FC6'    'FT10'    'T7'    'C3'...
            'CZ'    'C4'    'T8'    'CP5'    'CP1'    'CP2'    'CP6'    'P7'...
            'P3'    'PZ'    'P4'    'P8'    'O1'    'OZ'    'O2'};
switch chtype
    case 'Def'
        chlabels =  param.Ch32;
    case 'AR'
        chlabels = {'FP1'    'FPZ'    'FP2'    'F7'    'F3'    'FZ'    'F4'    'F8' ...
            'FC5'    'FC1'    'FC2'    'FC6'           'C3'...
            'CZ'    'C4'        'CP5'    'CP1'    'CP2'    'CP6'    ...
            'P3'    'PZ'    'P4'        'O1'    'OZ'    'O2'};
    case 'ARm'
        chlabels = { 'F7'    'F3'    'FZ'    'F4'    'F8' ...
            'FC5'    'FC1'    'FC2'    'FC6'           'C3'...
            'CZ'    'C4'        'CP5'    'CP1'    'CP2'    'CP6'    ...
            'P3'    'PZ'    'P4'        'O1'    'OZ'    'O2'};
end

param.ChInitial                 = chlabels;
param.NumChIni                  = length(param.ChInitial);
param.Ch                        = chlabels;
param.NumCh                     = length(param.Ch);

%% BCI param
param.repeat = 10;
param.Stims                     = 1 : 4;
param.NumStims                  = size(param.Stims,2);

%% Experiment param
param.NumTrTrial                  = 50;
param.NumTeTrial                  = 30;

param.switch_on                 = false;
param.Fs = 500;
%% epoching param
param.EpocType                  = 'singletr'; %{options:'meantr','meantr_tsb','singletr','singletr_tsb'}
param.Epoctime                  = 0.6;
param.Basetime                  = 0.2;
param.Epocline                  = param.Epoctime * param.Fs;
param.Baseline                  = param.Basetime * param.Fs;
param.Totalepoc                 = param.Epocline + param.Baseline;
param.Time                      = -param.Basetime : 1/param.Fs : param.Epoctime-1/param.Fs;
param.Numtrial                 = 0;

param.Sys                       = [12 13];
param.Targets                   = [];

%% figure
% param.H                         = figure(1); clf;
% set(param.H, 'color', 'w');
% 
% for i = 1:length(param.ChInitial)
%     param.SH(i)                       = axes;
%     hold(param.SH(i),'on');
% 
%     subplot(6,6,i,param.SH(i));
% 
%     param.h(i,1)                        = plot(nan,nan, 'parent',param.SH(i));
%     param.h(i,2)                        = plot(nan,nan, 'parent',param.SH(i));
%     param.h(i,3)                        = plot(nan,nan, 'parent',param.SH(i));
%     param.h(i,4)                        = plot(nan,nan, 'parent',param.SH(i));
%     param.h(i,5)                        = plot(nan,nan, 'parent',param.SH(i));
%     param.h(i,6)                        = plot(nan,nan, 'parent',param.SH(i));
% 
% end
% 
% set(gcf, 'Position', [200, 50, 1800, 1000])
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 10])

%% filter coeff
param.Fs                        = 500;
param.dFs                       = 125;
%-- IIR butterworth
[pB, pA]                        = butter(2, [0.5 1]./(param.Fs/2),'bandpass');
param.pBF                       = {pB, pA};

[B,A]                           = butter(4, [0.5 12]./(param.Fs/2),'bandpass');
param.BF                        = {B, A};

[lB,lA]                           = butter(4, [12]./(param.Fs/2),'low');
param.LF                        = {lB, lA};
[hB,hA]                           = butter(4, 0.5./(param.Fs/2),'high');
param.HF                        = {hB, hA};
[dB, dA]                       = butter(4, [50]./(param.Fs/2), 'low');
param.dLF                       = {dB, dA};

%-- FIR
B_fir = getfirfiltcoeff(0.5,'high',param.Fs,0);
param.HF_fir = {B_fir,1};

B_fir_lf = getfirfiltcoeff(50,'low',param.Fs,0);
param.dLF_fir = {B_fir_lf,1};

B_fir_lf_2 = getfirfiltcoeff(12,'low',param.Fs,0);
param.LF_fir = {B_fir_lf_2,1};

B_fir_ext = getfirfiltcoeff(0.5,'high',param.Fs,1);
param.HF_fir_ext = {B_fir_ext,1};

%% Preprocessing param
param.prep_factor = [1:4];

param.filterType = 'IIR';

param.DoAsr = 'Y';
param.windowlen = 0.5;
param.stepsize = floor(param.Fs*param.windowlen/2);
param.cutoff = 7;
param.maxdims = 0.66;

param.tmp.data = [];
param.tmp.srate = param.Fs;
param.tmp.cutoff = 1;

param.calibrate = false;
param.caltime = 60*param.Fs;

%% Classification param
param.FeatureType = 'DSP'; % {options: 'waveform','DSP'}
param.MovAvg = 0;          % {options: 1 (on), 0 (off)} 
param.trD.mode              = 'training'; % 'training' : training decoder , 'testing' : testing decoder

param.prediction                      = []; % predicted result
param.trD.kfold                     = 10;


end

function b = getfirfiltcoeff(cutoff,type,fs,extreme)
% editted EEGLAB pop_eegfiltnew
TRANSWIDTHRATIO = 0.25;

switch type
    case 'low'
        revfilt = 0;
    case 'high'
        revfilt = 1;
end
maxTBW = cutoff;

if revfilt == 0 % band/lowpass
    maxTBW(end) = fs/2 - cutoff(end);
elseif length(cutoff) == 2 % bandstop
    maxTBW = diff(cutoff)/2;
end
maxDf = min(maxTBW); % e.g. 0.5Hz higpass -> 0.5, 50Hz lowpass -> 50

if revfilt == 1
    df = min([max([maxDf*TRANSWIDTHRATIO 2]) maxDf]);
else
    df = min([max([cutoff(1) * TRANSWIDTHRATIO 2]) maxDf]);
end


filtorder = 3.3 / (df / fs); % Hamming window
filtorder= ceil(filtorder / 2) * 2; % Filter order must be even.

if extreme == 1
    filtorder = 1500;
end

winArray = windows('hamming', filtorder + 1);

switch type
    case 'low'
        cutoffArray = cutoff + df/2;
        b = firws(filtorder, cutoffArray / (fs/2), winArray);
        
    case 'high'
        cutoffArray = cutoff - df/2;
        b = firws(filtorder, cutoffArray / (fs/2), type, winArray);
        
end

end

