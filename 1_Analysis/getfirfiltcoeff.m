
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