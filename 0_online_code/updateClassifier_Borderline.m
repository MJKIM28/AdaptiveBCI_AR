function param = updateClassifier(Feature,output_init,param)

Nt = size(Feature,1)/(param.repeat*param.NumStims);
Nfeat = size(Feature,2);

%-- assign label to data. labels from previous classifier 
Feature_re = reshape(Feature',Nfeat,param.repeat,param.NumStims,Nt);
Feat_tar = Feature_re(:,:,output_init);
Feat_tar = Feat_tar(:,:)';
Feat_nar = Feature_re(:,:,setdiff(param.Stims,output_init));
Feat_nar = Feat_nar(:,:)';

Feature_new = [Feat_tar;Feat_nar];
label = [ones(size(Feat_tar,1),1); -ones(size(Feat_nar,1),1)];

%-- get score of current data
[sc] = predict_mj(param.trD,Feature_new);
sc_tar = sc(label==1,end);
sc_nar = sc(label==-1,end);

m = 10; k = 5;
%-- get danger set
feat_tar_old = param.trD.SV(param.trD.Y_SV == 1,:);
[eIdx,eD] = knnsearch(param.trD.SV,feat_tar_old ,'K', m+1, 'Distance', 'euclidean');
[eIdx,eD] = knnsearch(param.trD.SV,feat_tar_old ,'K', m+1, 'Distance', @distfun);

mysigmoid(param.trD.SV(end,:),param.trD.SV(end,:))
DANGER_set = find(sum(param.trD.Y_SV(eIdx(:,2:end))==-1, 2) > ((m+1)/2)); % find majority (over half) of nearest neighor is nontarget 
NInd = size(DANGER_set,1);

[Index, D] = knnsearch(Feature_new, feat_tar_old(DANGER_set,:),'K',k, 'Distance','euclidean');

samplesizenew = NInd*k ;

%-- get data on, within the margin, or missclassified ( t_n*y(x_n) <= 1 )
sv1 = find(sc_tar <= 1 ); 
Nsv = length(sv1);
sv2 = find(-sc_nar <= 1 );
Nsv2 = length(sv2);


Nsvall = Nsv + Nsv2;

%-- get alpha for new support vectors
alpha1 = []; alpha2 = [];
score_sv_tar = param.trD.score_SV(param.trD.Y_SV==1,end);
score_sv_nar = param.trD.score_SV(param.trD.Y_SV==-1,end);
alpha_tar = param.trD.Alpha(param.trD.Y_SV == 1);
alpha_nar = param.trD.Alpha(param.trD.Y_SV == -1);
for ii = 1:Nsv
    [~,id] = min(abs(score_sv_tar - sc_tar(sv1(ii))));
    alpha1(ii) = alpha_tar(id);
end
for ii = 1:Nsv2
    [~,id] = min(abs(score_sv_nar - sc_nar(sv2(ii))));
    alpha2(ii) = alpha_nar(id);
end
Alphatemp = [alpha1 alpha2]';
SV_Y =  [ones(Nsv,1);-ones(Nsv2,1)]; 


Alpha_new_temp = [param.trD.Alpha; Alphatemp];
labels = [param.trD.Y_SV;SV_Y];
SumAlpha = sum(Alpha_new_temp.*labels);
Alpha_new = Alpha_new_temp;

while abs(SumAlpha) > 0.0001 
    NotMaxAlpha1 = find(Alpha_new <1 & labels == 1);
    NotMaxAlpha2 = find(Alpha_new <1 & labels == -1);
    NNotMaxAlpha = length(NotMaxAlpha1) + length(NotMaxAlpha2);
    Alpha_new(NotMaxAlpha1) = Alpha_new(NotMaxAlpha1) - SumAlpha/NNotMaxAlpha;
    Alpha_new(NotMaxAlpha2) = Alpha_new(NotMaxAlpha2) + SumAlpha/NNotMaxAlpha;

    Alpha_new(Alpha_new < 0) = 0;
    Alpha_new(Alpha_new > 1) = 1;

    SumAlpha = sum(Alpha_new.*labels);
    fprintf('.');
end
fprintf('\n');

SV_new = [param.trD.SV; Feat_tar(sv1,:);Feat_nar(sv2,:)];
Y_SV_new = [param.trD.Y_SV; SV_Y];

zeroalpha = find(Alpha_new == 0);

Alpha_new(zeroalpha) = [];
SV_new(zeroalpha,:) = [];
Y_SV_new(zeroalpha) = [];

[~,ida] = min(param.trD.Alpha);
Bias_new = mean(Y_SV_new(Alpha_new < 1,:) - sum(mysigmoid(SV_new(Alpha_new > 0,:),SV_new(Alpha_new < 1,:)).*Y_SV_new(Alpha_new > 0,:).*Alpha_new(Alpha_new > 0,:))');

score_SV_new = [param.trD.score_SV; [-sc_tar(sv1,:) sc_tar(sv1,:)]; [sc_nar(sv2,:) -sc_nar(sv2,:)]];
score_SV_new(zeroalpha,:) = [];

param.trD.Alpha = Alpha_new;
param.trD.Biad =  Bias_new;
param.trD.SV = SV_new;
param.trD.Y_SV = Y_SV_new;
param.trD.score_SV = score_SV_new;
end