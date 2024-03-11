function [EP,param]     = Epoching_singletr_TSBAR(sig, trig, param)
try
    latency         = find(trig ~= 0); % Find latency of triggers
    type            = trig(trig ~= 0); % Find triggers
    
    init_t         = find(type == 11);
    start_t         = find(type == 12);
    end_t           = find(type == 13 );%%| type == 14);
    
    init_len = length(init_t);
    start_len = length(start_t);
    end_len = length(end_t);
    
    target_set = repmat([1 2 3 4]',20,1);
    target_use = target_set(1:end_len);
   
    
    
    target = zeros(end_len,1);
    START = zeros(end_len,1);
    
    %%
    %     if strcmp(param.decoder.mode,'training')
    %     ForBlockNum = [init_len start_len end_len];
    %     [~,w]= max(ForBlockNum);
    %     NumBlockEstimated = ForBlockNum(w);
    %
    %     % 11,12,13 중에 어느 거 기준으로  epoching할지 정하기: 11 고정
    % %     switch w
    % %         case 1
    % %             DivideBlock = init_t;
    % %         case 2
    % %             DivideBlock = start_t;
    % %         case 3
    % %             DivideBlock = end_t;
    % %     end
    %     DivideBlock = init_t;
    %     START = DivideBlock;
    %
    %     % 빠진  block 찾기
    %     Missed =   find(diff(DivideBlock) > 50) + 1;
    %     target_temp = target_set([1:Missed-1,Missed+1:end]);
    %     target_use = target_temp(1:NumBlockEstimated);
    %
    %     % epoching 하기
    %
    %
    %     elseif strcmp(param.decoder.mode,'testing')
    %         NumBlockEstimated = 1;
    %           START(1) = 1;
    % end
    %%
    if strcmp(param.trD.mode,'training')
        % find first target
        %         if type(1) == 10 %% trigger "10" in first trigger
        %             type(1) = []; latency(1) = [];
        %         end
        if type(1+1) == 11 && type(3+1) == 12 % 1st block, 11 and 12
%         if type(1) == 11 && type(2) == 12 % 1st block, 11 and 12
            target(11) = type(2+1);
            START(1) = 3+1;
        elseif type(2+1) == 12 % 1st block, only 12 (without 11)
            target(1) = type(1+1);
            START(1) = 2+1;
        elseif type(1+1) == 11 && latency(2+1) - latency(1+1) < 0.2*param.Fs % 1st block, only 11 (without 12)
            target(1) = type(2+1);
            START(1) = 2+1;
        elseif ismember(type(1+1),param.Stims) && ismember(type(2+1),param.Stims) % 1st block, both 11 and 12 ignored
            if latency(2+1) - latency(1+1) < 0.1*param.Fs
                target(1) = type(1+1);
                START(1) = 1+1;
            end
        end
        
        
        if (init_len ~= param.NumTrTrial) || (start_len ~= param.NumTrTrial) % at least one of 11s or 12s ignored
            for t = 1:end_len-1
                if type(end_t(t)+1) == 11 % if init (11) detected in t+1 th block
                    if type(end_t(t)+3) == 12 % if the one after the next of init (11) is start (12) (e.g. 13 11 2 12)
                        target(t+1) = type(end_t(t)+2);
                        START(t+1) = end_t(t)+3;
                    elseif latency(end_t(t)+2) - latency(end_t(t)+1) < 0.02*param.Fs % start (12) omitted (e.g. 13 11 2 3)
                        target(t+1) = type(end_t(t)+2);
                        START(t+1) = end_t(t)+2;
                    end
                elseif ismember(type(end_t(t)+1),param.Stims) % (1~4) right after end(13)
                    if type(end_t(t)+2) == 12 % if start(12) is after (1~4) (e.g. 13 2 12 3)
                        target(t+1) = type(end_t(t)+1);
                        START(t+1) = end_t(t)+2;
                    elseif latency(end_t(t)+2) - latency(end_t(t)+1) < 0.1*param.Fs %(1~4) after (1~4), interval is less than 100ms (e.g. 13 2 3)
                        target(t+1) = type(end_t(t)+1);
                        START(t+1) = end_t(t)+1;
                    end
                end
            end
        else
            START = start_t;
        end
        singleTr = cell(param.NumStims,param.repeat);
    elseif strcmp(param.trD.mode,'testing')
        if init_len ~= 1 || start_len ~= 1
            START(1) = 1;
        else
            START = start_t;
        end
    end
    
    %%
    
%     target_use = [target_use;  2];
%     target_use([33]) = [];
%     START([32]) = []; %end_t([25 26])= [];

    end_len = length(end_t);
    
    
    numTr           = end_len;
    %     numTr = NumBlockEstimated;
    
    tar                 = NaN(param.NumCh, param.Totalepoc, numTr);
    nar                 = NaN(param.NumCh, param.Totalepoc, 3, numTr);
    if ~strcmp(param.trD.mode,'training')
        nar                 = zeros(param.NumCh, param.Totalepoc, 4, numTr);
    end
    
%     tmp = NaN(param.NumCh, param.Totalepoc,param.NumRepeat,param.NumStims,numTr);
EP.dat = NaN(param.Totalepoc,param.NumCh,param.NumStims*param.repeat,numTr);
EP.lat = NaN(param.Totalepoc,param.NumStims*param.repeat,numTr);

for b = 1 : numTr
        if START(b) ~= 0
            %         if b < numTr
            latency_onset            = latency(START(b)+1 : end_t(b)-1);
            Diff = diff(latency_onset)/param.Fs;
            SumDiff = movsum(Diff,2);
            DelayedIdx = find(Diff < 0.15);
            IsDelayed = SumDiff(DelayedIdx) < 0.6;
            latency_onset(DelayedIdx(IsDelayed)) = NaN;
            
            
%              (latency(end_t(b)) - latency(START(b)))./500
%             stim_type           = type(START(b)+1 : end_t(b)-1);
%             stim_type(DelayedIdx(IsDelayed)) = NaN;
            
             
            
    
            baseline = latency_onset - param.Baseline;
            epoc = latency_onset + param.Epocline - 1;
            
            for n = 1:length(baseline)
                if ~isnan(baseline(n)) && ~isnan(epoc(n))
                    EP.dat(:,:,n,b) = sig(:,baseline(n):epoc(n))';
                    EP.lat(:,n,b) = trig(baseline(n):epoc(n));
                end
            end
    
    
%             for s = 1 : param.NumStims
%                 ind             = stim_type == param.Stims(s);
%                 trials          = latency_onset(ind);
%                 for t = 1 :length(trials)
%                     ind_sig     = trials(t) - param.Baseline : trials(t) + param.Epocline - 1;
%                     tmp(:,:,t,s,b) = sig(:,ind_sig);
%                     
%                 end
%             end
        end
    end
    
    
%%%%% remove bad trial
%     if strcmp(param.decoder.mode,'training')
%     [tmp,rmtr, param.badtrth] = removeoutlier(tmp,500,[]);
%     else
%         [tmp,rmtr, param.badtrth] = removeoutlier(tmp,500,param.badtrth);    
%     end
%     figure(20); plot(squeeze(tmp(1,:,:)))
%%%%%%%%%
    
% erp_tmp         = squeeze(sum(tmp,3,'omitnan')./ length(trials));
%     erp      = (erp_tmp - repmat(mean(erp_tmp(:, 1 : param.Baseline,:,:),2 ,'omitnan'), 1, size(erp_tmp,2), 1, 1))./...
%         repmat(std(erp_tmp(:, 1 : param.Baseline,:,:),[],2,'omitnan'), 1, size(erp_tmp,2),1,1);
%     
%     if ~strcmp(param.decoder.mode,'training')
%         nar(:,:,:,b)    = erp;
%     else
%         for b = 1:numTr
%             tar(:,:,b)      = erp(:,:,(target_use(b) == param.Stims),b);
%             nar(:,:,:,b)    = erp(:,:,~(target_use(b) == param.Stims),b);
%             
%         end
%     end
%     
%     
%     ch = find(param.ChDraw);
%     for k = 1:param.NumCh
%         if strcmp(param.decoder.mode,'training')
%         set(param.h(k,1), 'XData', param.Time,'YData', squeeze(mean(tar(k,:,:),3,'omitnan')), 'Color','r','LineWidth',2);
%         set(param.h(k,2), 'XData', param.Time,'YData', squeeze(mean(nar(k,:,1,:),4,'omitnan')),'Color','g');
%         set(param.h(k,3), 'XData', param.Time,'YData', squeeze(mean(nar(k,:,2,:),4,'omitnan')),'Color','b');
%         set(param.h(k,4), 'XData', param.Time,'YData', squeeze(mean(nar(k,:,3,:),4,'omitnan')),'Color','k');
%         else 
%             set(param.h(k,1), 'XData', param.Time,'YData', squeeze(mean(nar(k,:,1,:),4,'omitnan')),'Color','r','LineWidth',1);
%             set(param.h(k,2), 'XData', param.Time,'YData', squeeze(mean(nar(k,:,2,:),4,'omitnan')),'Color','g','LineWidth',1);
%             set(param.h(k,3), 'XData', param.Time,'YData', squeeze(mean(nar(k,:,3,:),4,'omitnan')),'Color','b','LineWidth',1);
%             set(param.h(k,4), 'XData', param.Time,'YData', squeeze(mean(nar(k,:,4,:),4,'omitnan')),'Color','k','LineWidth',1);
%         end
%     end
     
%     EP.tar = tar;
%     EP.nar = nar;
    EP.target = target_use;
    %     EP.single = singleTr;
    
    clearvars -except EP param
catch
    keyboard
end