function [Feature,labels, param] = FeatureExt_basic(EP, param)
% try
NumBlocks                   = size(EP.tar,3);
Win = 1;
ind = param.Baseline+0.15*param.Fs+1:Win:param.Totalepoc;
if param.MovAvg 
    MovAvgWin = 10;
MovAvgDS = fix(MovAvgWin/2);
else 
    MovAvgWin = 1;
    MovAvgDS = 1;
end
if strcmp(param.trD.mode,'training')
        tar_temp                = permute(EP.tar(:,ind,:),[2,1,3]);
        tar_temp = movmean(tar_temp,MovAvgWin,1);
        tar_temp = tar_temp(1:MovAvgDS:end,:,:);
        tar_temp2                = reshape(tar_temp,size(tar_temp,1)*size(tar_temp,2),size(tar_temp,3));
        tar_feat = tar_temp2';
        
        nar_temp = permute(EP.nar(:,ind,:,:),[2,1,3,4]);
        nar_temp = movmean(nar_temp,MovAvgWin,1);
        nar_temp = nar_temp(1:MovAvgDS:end,:,:,:);
        nar_temp2 = reshape(nar_temp,size(nar_temp,1)*size(nar_temp,2),size(nar_temp,3)*size(nar_temp,4));
        nar_feat= nar_temp2';
        
        labels                = [ones(size(tar_feat,1),1); -ones(size(nar_feat,1),1)];
        
        Feature = [tar_feat; nar_feat];
  
else
        feat_temp = permute(EP.nar(:,ind,:,:),[2,1,3,4]);
        feat_temp = movmean(feat_temp,MovAvgWin,1);
        feat_temp = feat_temp(1:MovAvgDS:end,:,:,:);
        feat_temp2 = reshape(feat_temp,size(feat_temp,1)*size(feat_temp,2),size(feat_temp,3)*size(feat_temp,4));
        Feature= feat_temp2';
        labels = [];
 

  
end
clearvars -except Feature labels param