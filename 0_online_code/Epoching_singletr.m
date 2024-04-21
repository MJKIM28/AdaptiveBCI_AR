function EP = Epoching_singletr(signal,trigger,param)
% EP.dat: single epochs 
%       ( epoch_length x channel_size x (stim_type_size x repetition_number) x block_size ) 
% EP.lat: event markers in same size of epochs
%       ( epoch_length x (stim_type_size x repetition_number) x block_size ) 

try

latency = find(trigger ~= 0);
type = trigger(trigger ~= 0);

ind_start = find(type == param.Sys(1)); % trigger 12: start of a block
ind_end = find(type == param.Sys(2));   % trigger 13: end of a block

% Re-check the start and the end of a block
if length(ind_start) ~= length(ind_end)
    % Between the start and the end of a block,
    % # of trigger should be stim_type_size x repetition_number + 1 
   getdiff = bsxfun(@minus,ind_end,ind_start');  
   [r,c] = find(getdiff == param.NumStims*param.repeat+1); 
    
   ind_start = ind_start(r); 
   ind_end = ind_end(c);   
end

Nb = length(ind_start);
EP.target = type(ind_start - 1); 

EP.dat = NaN(param.Totalepoc,param.NumCh,param.NumStims*param.repeat,Nb);
EP.lat = NaN(param.Totalepoc,param.NumStims*param.repeat,Nb);
for b = 1:Nb
    ind_block = ind_start(b)+1:ind_end(b)-1;
    latency_onset = latency(ind_block);
    
    baseline = latency_onset - param.Baseline;
    epoc = latency_onset + param.Epocline - 1;
    
    for n = 1:length(baseline)
        EP.dat(:,:,n,b) = signal(:,baseline(n):epoc(n))';
        EP.lat(:,n,b) = trigger(baseline(n):epoc(n));
    end
    
end

catch
    keyboard
end