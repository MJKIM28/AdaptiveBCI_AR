function    ouput_st = getoutputStatic(EP_1block,param,useUF,useLinear)
if ~useUF % if not use updated feature
    param.DSP = param.DSP.init;
end

[Feature,label,param] = FeatureExtraction_newCV(EP_1block,param);

global gamma
gamma = param.trD.gamma;

if ~useLinear
    [C,sc] = predict(param.trD.mdl_init, Feature);

    Nt = size(Feature,1)/(param.repeat*param.NumStims);
    tescore_bl = reshape(sc(:,end),param.repeat,param.NumStims,Nt);
    tescore_bl_sum = sum(tescore_bl,1);
    [~,hit_id] = max(tescore_bl_sum,[],2);
    ouput_st = squeeze(hit_id);
else
    Nt = size(Feature,1)/(param.repeat*param.NumStims);
    Feat_re = reshape(Feature',size(Feature,2),param.repeat,param.NumStims,Nt);

    Feat_use =  reshape(mean(Feat_re,2),size(Feature,2),param.NumStims,Nt)'; % Nt = 1
    [C,sc] = predict(param.trD.mdl_init, Feat_use);
    [~,hit_id] = max(sc(:,end));
    ouput_st = squeeze(hit_id);

end

end