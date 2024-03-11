
%% ICL_feature_extractor
function features = ICL_feature_extractor_MJ(act, winv, chanloc, param)
ncomp = size(act, 1);
pnts = size(act, 2);
%% calc topo
topo = zeros(32, 32, 1, ncomp);
for it = 1:ncomp
        [~, temp_topo, plotrad] = ...
            topoplotFast(winv(:, it), chanloc, ...
            'noplot', 'on');
  
    temp_topo(isnan(temp_topo)) = 0;
    topo(:, :, 1, it) = temp_topo / max(abs(temp_topo(:)));
end

% cast
topo = single(topo);
    
%% calc psd
psd = eeg_rpsd_MJ(act, 100, param);

% extrapolate or prune as needed
nfreq = size(psd, 2);
if nfreq < 100
    psd = [psd, repmat(psd(:, end), 1, 100 - nfreq)];
end

% undo notch filter
for linenoise_ind = [50, 60]
    linenoise_around = [linenoise_ind - 1, linenoise_ind + 1];
    difference = bsxfun(@minus, psd(:, linenoise_around), ...
        psd(:, linenoise_ind));
    notch_ind = all(difference > 5, 2);
    if any(notch_ind)
        psd(notch_ind, linenoise_ind) = mean(psd(notch_ind, linenoise_around), 2);
    end
end

% normalize
psd = bsxfun(@rdivide, psd, max(abs(psd), [], 2));

% reshape and cast
psd = single(permute(psd, [3 2 4 1]));

%% autocorr

if pnts / param.Fs > 5
            autocorr = eeg_autocorr_welch_MJ(act, param);
        else
            autocorr = eeg_autocorr_MJ(act, param);
        end
    autocorr = single(permute(autocorr, [3 2 4 1]));

features = {0.99 * topo, 0.99 * psd, 0.99 * autocorr};
