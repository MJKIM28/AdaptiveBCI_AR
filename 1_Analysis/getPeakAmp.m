function [Amplitude,PeakLat] = getPeakAmp(ERPTotal,Win,winsize)
% ERPTotal: ch x Time
% Amplitude: ch x (win type) x (negative/positive)
Nch = size(ERPTotal,1);
for w = 1:length(Win)
    for ch = 1:Nch
        [peakval_pos,peakind_pos] = findpeaks(ERPTotal(ch,Win{w})); % positive peaks
        [peakval_neg,peakind_neg] = findpeaks(-ERPTotal(ch,Win{w})); % negative peaks
        [maxval,maxind] = max(abs(ERPTotal(ch,Win{w}))); % maximum amplitude (absolute) within the window


        %-- negative peak
        peakval_neg_ = ERPTotal(ch,Win{w}(peakind_neg));
        [~,peakselid] = min(peakval_neg_);
        peaksellat = Win{w}(peakind_neg(peakselid));

        if ~isempty(peaksellat)
        PeakLat(ch,w,1) = peaksellat;
        WinUse = PeakLat(ch,w,1) - winsize: PeakLat(ch,w,1) + winsize;

        Amplitude(ch,w,1) = mean(ERPTotal(ch,WinUse),2,'omitnan');
        else
             PeakLat(ch,w,1) = NaN;
             Amplitude(ch,w,1) = NaN;
        end
             

        %-- positive peak
        [~,peakselid] = max(peakval_pos);
        peaksellat = Win{w}(peakind_pos(peakselid));
       if ~isempty(peaksellat)
        PeakLat(ch,w,2) = peaksellat;
        WinUse = PeakLat(ch,w,2) - winsize: PeakLat(ch,w,2) + winsize;

        Amplitude(ch,w,2) = mean(ERPTotal(ch,WinUse),2,'omitnan');
        else
             PeakLat(ch,w,2) = NaN;
             Amplitude(ch,w,2) = NaN;
        end


        %         if ismember(maxind,peakind_neg) || ismember(maxind,peakind_pos) % maximum value in positive or negative peaks
        %
        %             %               [~,peakselectid] = max(peakval);
        %             peaksellat = fix(Win{w}(maxind));
        %
        %         else % find maximum peak (negative or positive)
        %             [~,choose] = max([max(abs(peakval_neg)) max(abs(peakval_pos))]);
        %             if choose == 1
        %                 [~,peakselid] = max(abs(peakval_neg));
        %                 peaksellat = Win{w}(peakind_neg(peakselid));
        %
        %             elseif choose == 2
        %                 [~,peakselid] = max(abs(peakval_pos));
        %                 peaksellat = Win{w}(peakind_pos(peakselid));
        %
        %             end
        %         end
        %
        %         PeakLat(ch,w) = peaksellat;
        %         WinUse = PeakLat(ch,w) - winsize: PeakLat(ch,w) + winsize;
        %
        %         Amplitude(ch,:,:,w) = mean(ERPTotal(ch,WinUse,:,:),2,'omitnan');

    end
end
end

