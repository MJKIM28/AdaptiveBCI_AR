function param = updateClassifier_old(Feature,output_init,param)

Nt = size(Feature,1)/(param.repeat*param.NumStims);
Nfeat = size(Feature,2);

Feature_re = reshape(Feature',Nfeat,param.repeat,param.NumStims,Nt);
Feat_tar = Feature_re(:,:,output_init);
Feat_tar = Feat_tar(:,:)';
Feat_nar = Feature_re(:,:,setdiff(param.Stims,output_init));
Feat_nar = Feat_nar(:,:)';

Feature_new = [Feat_tar;Feat_nar];
label = [ones(size(Feat_tar,1),1); -ones(size(Feat_nar,1),1)];

[sc] = predict_mj(param.trD,Feature_new);
sc_tar = sc(label==1,end);
sc_nar = sc(label==-1,1);
svtemp1 = find(abs(sc_tar) <1);
Nsv = length(svtemp1);
[scsort,scid] = sort(sc_nar,'ascend');
svtemp2 = scid(1:Nsv);
sv2 = svtemp2(abs(sc_nar(svtemp2)) <= 1);
Nsv2 = length(sv2);

[scsort2,scid2] = sort(sc_tar(svtemp1),'ascend');
sv1 = svtemp1(scid2(1:Nsv2));

alpha1 = []; alpha2 = [];
score_sv_tar = param.trD.score_SV(param.trD.Y_SV==1,end);
score_sv_nar = param.trD.score_SV(param.trD.Y_SV==-1,end);
alpha_tar = param.trD.Alpha(param.trD.Y_SV == 1);
alpha_nar = param.trD.Alpha(param.trD.Y_SV == -1);
for ii = 1:Nsv2
    [~,id] = min(abs(score_sv_tar - sc_tar(sv1(ii))));
    alpha1(ii) = alpha_tar(id);

    [~,id] = min(abs(score_sv_nar - sc_nar(sv2(ii))));
    alpha2(ii) = alpha_nar(id);
end
Alphatemp = [alpha1 alpha2]';
SV_Y =  [ones(Nsv2,1);-ones(Nsv2,1)];

SumAlpha = sum(Alphatemp.*SV_Y);
if SumAlpha ~= 0
    NotMaxAlpha = find(Alphatemp <1);
    Alphatemp(NotMaxAlpha) = Alphatemp(NotMaxAlpha) - SumAlpha/length(NotMaxAlpha);
end
Alpha_new = [param.trD.Alpha; Alphatemp];

SV_new = [param.trD.SV; Feat_tar(sv1,:);Feat_nar(sv2,:)];
Y_SV_new = [param.trD.Y_SV; SV_Y];

[~,ida] = min(param.trD.Alpha);
Bias_new = mean(Y_SV_new(Alpha_new < 1,:) - sum(mysigmoid(SV_new(Alpha_new > 0,:),SV_new(Alpha_new < 1,:)).*Y_SV_new(Alpha_new > 0,:).*Alpha_new(Alpha_new > 0,:))');

score_SV_new = [param.trD.score_SV; [-sc_tar(sv1,:) sc_tar(sv1,:)]; [sc_nar(sv2,:) -sc_nar(sv2,:)]];
param.trD.Nadded = length(sv1) + length(sv2);
param.trD.Alpha = Alpha_new;
param.trD.Bias =  Bias_new;
param.trD.SV = SV_new;
param.trD.Y_SV = Y_SV_new;
param.trD.score_SV = score_SV_new;
param.trD.score_new = {sc_tar,-sc_nar};

 
end