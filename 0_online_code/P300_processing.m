function [C, param]   = P300_processing(sig, trig, param)
global DS
% updated @ 20210426
try
if isfield(param.trD, 'mdl')
    param.trD.mode  = 'testing';
    fprintf('Testing start!..\n');
else
    fprintf('Training start!..\n');
end

sig = double(sig*0.04883);
[sig, param]            = PreProcess(sig,param);

EP = EpochData(sig,trig,param);
[Feature, label, param] = FeatureExtraction(EP, param);



[C,param]               = Classification(Feature,label, param);

% FOR Dynamic Stopping
% If iteration end but block not end,
% decide to stop or not
if (~param.switch_on) && isfield(param,'switch_on_iter')
   dat = reshape(param.decoder.data,param.NumStims,size(param.decoder.data,1)/param.NumStims);
   
   DS.data = dat;
   DS = getpval(DS);
   [param.DS.stop,param.DS.class] = decidestopping(DS);
end

% clearvars -except C param
catch 
    keyboard
end