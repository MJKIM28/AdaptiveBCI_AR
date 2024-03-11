function psdmed = eeg_rpsd_MJ(act, nfreqs, param, pct_data)

% clean input cutoff freq
nyquist = floor(param.Fs / 2);

if nfreqs > nyquist
    nfreqs = nyquist;
end
if ~exist('pct_data', 'var') || isempty(pct_data)
    pct_data = 100;
end

% setup constants
ncomp = size(act, 1);
pnts = size(act,2);
n_points = min(size(act,2), param.Fs);
window = windows('hamming', n_points, 0.54)';
cutoff = floor(pnts / n_points) * n_points;
index = bsxfun(@plus, ceil(0:n_points / 2:cutoff - n_points), (1:n_points)');

if ~exist('OCTAVE_VERSION', 'builtin')
    rng('default')
else
    rand('state', 0);
end

n_seg = size(index, 2) ;
subset = randperm(n_seg, ceil(n_seg * pct_data / 100)); % need to improve this

if exist('OCTAVE_VERSION', 'builtin') == 0
    rng('shuffle')
end

% calculate windowed spectrums
psdmed = zeros(ncomp, nfreqs);
for it = 1:ncomp
    temp = reshape(act(it, index, :), [1 size(index) .* [1 1]]);
    temp = bsxfun(@times, temp(:, :, subset), window);
    temp = fft(temp, n_points, 2);
    temp = temp .* conj(temp);
    temp = temp(:, 2:nfreqs + 1, :) * 2 / (param.Fs*sum(window.^2));
    if nfreqs == nyquist
        temp(:, end, :) = temp(:, end, :) / 2; end

    psdmed(it, :) = 20 * log10(median(temp, 3));
end
end
