function[hitLabel,param]= Classification_new(Feature,label,param)
gamma = param.trD.gamma;
%         gamma = 1/(size(Feature,2)*var(Feature(:)));
[sc] = predict_mj(param.trD, Feature);

Nt = size(Feature,1)/(param.repeat*param.NumStims);
tescore_bl = reshape(sc(:,end),param.repeat,param.NumStims,Nt);
tescore_bl_sum = sum(tescore_bl,1);
[~,hit_id] = max(tescore_bl_sum,[],2);
hit_id = squeeze(hit_id);
hitLabel = param.Stims(hit_id);

param.trD.score = tescore_bl;

end
