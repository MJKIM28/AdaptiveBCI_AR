function [sig,param] = PreProcess(sig,param)

%% 0.5Hz HPF 
if isfield(param,'filterType') && strcmp(param.filterType,'FIR')
    if strcmp(param.trD.mode,'training')
        b = param.HF_fir{1};
    elseif strcmp(param.trD.mode,'testing')
        b = param.HF_fir_ext{1};
    end
        sig = filtfilt(b,1,double(sig)')';   
else
    sig = filtfilt(param.HF{1}, param.HF{2}, double(sig)')'; %180213 updated: default filtering 0.5 Hz
end
fprintf('[1] 0.5Hz HPF .. Done \n')
param.NumCh = size(sig,1);

%% Bad channel rejection
% clf(figure(10));
if ismember(1,param.prep_factor) % 1: bad channel rejection
    if ~isfield(param,'badch') && ~strcmp(param.trD.mode,'testing')
        badch = prebadchannelrejection(sig,param);
        fprintf('Bad channels: %d\n', badch);
        param.badch = badch;
        for i = 1:length(param.badch)
            param.Ch{param.badch(i)} =[];
        end
        param.Ch = param.Ch(~cellfun('isempty',param.Ch));
    end
    
%     if ~isempty(param.badch)
%         figure(10); subplot(2,2,1); plot(sig(param.badch(1),:)); drawnow
%     else
%         figure(10); subplot(2,2,1); plot(sig(1,:)); drawnow
%     end
%     
    fprintf('[2] Bad channel Rejection .. Done \n');
  
    if isfield(param,'interp')
        if param.interp == 1
            sig = eeg_interp_MJ(sig,param.badch,'spherical');
            fprintf('[2]_1 Interpolation .. Done \n');
        else
            sig(param.badch,:) = [];
            param.NumCh = size(sig,1);
        end
    else
%         if ~isempty(param.badch)
%             hold on;
%             plot(sig(param.badch(1),:));
%         else
%             hold on;
%             plot(sig(1,:));
%         end
        sig(param.badch,:) = [];
         param.NumCh = size(sig,1);

    end
end

%% Rereferencing
if ismember(2,param.prep_factor) % 3: CAR & REST
%     figure(10); subplot(2,2,2); plot(sig(1,:)); hold on;drawnow
    sig = sig - repmat(mean(sig,1),size(sig,1),1);
    %     figure(10); subplot(2,2,3); plot(sig(1,:));drawnow
    fprintf('[3] CAR .. Done  \n');
    
    %%% REST %%%
    %     if isfield(param,'REST')
    %         if ~strcmp(param.decoder.mode,'testing') && ~isfield(param.REST,'G')
    %             if param.NumCh ~= 31  % bad channel rejected -> calculate leadfield matrix again
    %                 file = 'D:\[1]EEGBCI\[2]Research\Code&Algorithm\32_ch_XYZ_loc_Nose_Y.txt';
    %                 filenew = editchanloc(file,param.badch);
    %                 param.REST.calculateG = 1;
    %                 Gfile = [];
    %             else
    %                 Gfile =  'D:\[1]EEGBCI\[2]Research\Code&Algorithm\REST\Lead_Field_32ch.dat';
    %                 param.REST.calculateG = 0;
    %             end
    %             [sig,G] = REST(sig,[],Gfile,param.NumCh,param.REST.calculateG);
    %             param.REST.LeadField = G;
    %             param.REST.calculateG = 0;
    %                     fprintf('[3]_1 REST .. Done  \n');
    %
    %         else
    %
    %             [sig,~] = REST(sig,param.REST.LeadField,param.NumCh,param.REST.calculateG);
    %                     fprintf('[3]_1 REST .. Done  \n');
    %
    %         end
    %
    %     end
    hold on;
%     plot(sig(1,:));
    
    
end

%% 50Hz LPF : avoid line noise & high frequency noise
if isfield(param,'filterType') && strcmp(param.filterType,'FIR')
    %     figure;
    %     plot(sig(1,:))
    b = param.dLF_fir{1};
    sig = filtfilt(b,1,sig')';
    %     hold on; plot(sig(1,:));
else
    sig = filtfilt(param.dLF{1}, param.dLF{2}, double(sig)')'; %180213 updated: default filtering  50Hz
end
fprintf('[4] 50Hz LPF .. Done \n')

%% ASR
if ismember(3, param.prep_factor) % 3: ASR
%     figure(10); subplot(2,2,3); plot(sig(1,:)); hold on;drawnow
    if (isfield(param,'filterType') && strcmp(param.filterType,'FIR'))...
            ||(isfield(param,'cal') && strcmp(param.cal,'on'))
        % Not using 'cal_sig', but use ALL signal to get parameters
        ref = clean_windows_MJ(sig,param.Fs);
        param.state = asr_calibrate(ref,param.Fs,param.cutoff); 
        sigout = ASR(sig,param);
               
    else
        if ~strcmp(param.trD.mode,'testing')
%             load([param.dir,'/cal_sig.mat']);               
%             if size(cal_sig,1) ~= param.NumChIni
%                 ch_here = find(ismember(param.Ch32,param.ChInitial));
%                 cal_sig = cal_sig(ch_here,:); % follow param setting
%             end
%             param = preASR(cal_sig,param);
        ref = clean_windows_MJ(sig(:,1:5*param.Fs*60),param.Fs);
        param.state = asr_calibrate(ref,param.Fs,param.cutoff);
        
        end
        try
            switch param.DoAsr
                case 'Y'                 
                    sigout = ASR(sig,param);
                    fprintf('ASR is done..!\n');
                case 'N'
                    sigout = sig;
                    fprintf('Pass ASR ..!\n');
            end
        catch
            sigout = sig;
            fprintf('ASR is failed..!\n');
        end
    end
    
    sig = sigout;
%     figure(10); subplot(2,2,3); plot(sigout(1,:)); drawnow
    fprintf('[5] Artifact rejection .. Done \n')
    
end
%%
if ismember(5,param.prep_factor) % 5: ICA
    chanloc = readlocs('BP_channelLocs.locs');
    
    
    if length(chanloc) ~= param.NumChIni
        chrem = ~ismember({chanloc.labels},param.ChInitial);
        chanloc(chrem) = [];
    end
    if ~param.interp
        chanloc(param.badch) = [];
    end
    %%
    rank = getrank(sig,0);
    if isfield(param,'weights')
        w = param.weights*param.sphere;
    else
   
    [weights, sphere] = runica(sig,'pca',rank);
    w = weights*sphere;
    winv = pinv(w);

%      [winv, w] = fastica( sig, 'displayMode', 'off' ,'lastEig',rank);
%      weights = w;
%      sphere = eye(size(w,2));%No icasphere-> use eye() as icasphere in pop_runica 


%      winv    = sobi( sig);
%      weights = pinv(winv);
%      sphere = eye(size(weights,2));
%      w = weights*sphere;

    end
    
    act = w*sig;
    ncomp = size(act,1);
    features = ICL_feature_extractor_MJ(act,winv, chanloc, param);
    iclabels = run_ICL_MJ('default', features{:});
    
    [prop,lb] = max(iclabels,[],2);
    temp_rej_eye = find(lb == 3); % all eyes
    temp_rej = find(prop > 0.75);
    rem = find(lb(temp_rej) == 2 | lb(temp_rej) == 3 | ...
        lb(temp_rej) == 4 | lb(temp_rej) == 5 | lb(temp_rej) == 6); % 2: muscle, 3: eye, 4: heart , 5: line noise, 6: channel noise
    ic_rej = unique([temp_rej(rem);temp_rej_eye]);
    
    ic_use = setdiff(1:ncomp,ic_rej);
    sig = winv(:,ic_use)*act(ic_use,:);
    param.ica.weights = weights;
    param.ica.sphere = sphere;
    param.ica.rmcomp = ic_rej;
    param.ica.rmcompinfo = [prop(ic_rej) lb(ic_rej)];
end  


%% 12Hz LPF
if ismember(4,param.prep_factor) % 4: filtering
%     figure(10); subplot(2,2,4); plot(sig(1,:)); hold on;drawnow
    
    if isfield(param,'filterType') && strcmp(param.filterType,'FIR')
        b = param.LF_fir{1};
        sig = filtfilt(b,1,sig')';
    else
        sig = filtfilt(param.LF{1},param.LF{2},sig')';
    end
%     figure(10); subplot(2,2,4); plot(sig(1,:)); drawnow
    fprintf('[6] 12Hz LPF .. Done \n');
end

param.ChDraw = 1:param.NumCh;


end

%% from eeglab
function tmprank2 = getrank(tmpdata, pca_opt)
    
    tmprank = rank(tmpdata);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Here: alternate computation of the rank by Sven Hoffman
    %tmprank = rank(tmpdata(:,1:min(3000, size(tmpdata,2)))); old code
    covarianceMatrix = cov(tmpdata', 1);
    [~, D] = eig (covarianceMatrix);
    rankTolerance = 1e-7;
    tmprank2=sum (diag (D) > rankTolerance);
    if tmprank ~= tmprank2
        if pca_opt ~= 0
            fprintf('Warning: fixing rank computation inconsistency (%d vs %d) most likely because running under Linux 64-bit Matlab\n', tmprank, tmprank2);
        end
        tmprank2 = max(tmprank, tmprank2);
    end
end

