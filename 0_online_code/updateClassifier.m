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

if strcmp(param.trD.adaptmode,'margin')
    %-- get data on, within the margin, or missclassified ( t_n*y(x_n) <= 1 )
    sv1 = find(sc_tar <= 1.09 );
    Nsv = length(sv1);
    sv2 = find(-sc_nar <= 1.09 );
    Nsv2 = length(sv2);

    Nsvall = Nsv + Nsv2;

    %-- get alpha for new support vectors
    alpha1 = []; alpha2 = [];
    NoneOneAlpha = find(sc_tar(sv1) > 0.9 & sc_tar(sv1) < 1.1);
    alpha1(NoneOneAlpha) = 0;
    alpha1(setdiff(1:Nsv,NoneOneAlpha)) = 1;

    NoneOneAlpha = find(sc_nar(sv2) > 0.9 & sc_nar(sv2) < 1.1);
    alpha2(NoneOneAlpha) = 0;
    alpha2(setdiff(1:Nsv2,NoneOneAlpha)) = 1;

else
    %-- all samples
    sv1 = 1:size(sc_tar,1);
    Nsv = length(sv1);
    sv2 = 1:size(sc_nar,1);
    Nsv2 = length(sv2);

    Nsvall = Nsv + Nsv2;

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

end
Alphatemp = [alpha1';alpha2'];
SV_Y = [ones(length(alpha1),1);-ones(length(alpha2),1)];
SV_add = [Feat_tar(sv1,:);Feat_nar(sv2,:)];
sc_add = [[-sc_tar(sv1,:) sc_tar(sv1,:)]; [-sc_nar(sv2,:) sc_nar(sv2,:)]];

%-- set balanced data 
Ntar = length(find(param.trD.Y_SV == 1));
Nnar = length(find(param.trD.Y_SV == -1));

if Ntar < Nnar
    %-- only add tar SV
    Nsvnew = min([Nnar - Ntar Nsv]);
    Alphatemp = alpha1(1:Nsvnew)';
    SV_Y = ones(Nsvnew,1);
    SV_add = Feat_tar(sv1(1:Nsvnew),:);
    sc_add = [-sc_tar(sv1(1:Nsvnew),:) sc_tar(sv1(1:Nsvnew),:)];
elseif Ntar > Nnar
%     keyboard
    %-- add same number
    DelTar = Ntar - Nnar;
    Nsvnew = min([Nsv Nsv2]);
    Alphatemp = [alpha1(1:Nsvnew-DelTar) alpha2(1:Nsvnew)]';
    SV_Y =  [ones(Nsvnew-DelTar,1);-ones(Nsvnew,1)];
    SV_add = [Feat_tar(sv1(1:Nsvnew-DelTar),:);Feat_nar(sv2(1:Nsvnew),:)];
    sc_add = [[-sc_tar(sv1(1:Nsvnew-DelTar),:) sc_tar(sv1(1:Nsvnew-DelTar),:)]; [-sc_nar(sv2(1:Nsvnew),:) sc_nar(sv2(1:Nsvnew),:)]];

    


else
    %-- add same number
    Nsvnew = min([Nsv Nsv2]);
    Alphatemp = [alpha1(1:Nsvnew) alpha2(1:Nsvnew)]';
    SV_Y =  [ones(Nsvnew,1);-ones(Nsvnew,1)];
    SV_add = [Feat_tar(sv1(1:Nsvnew),:);Feat_nar(sv2(1:Nsvnew),:)];
    sc_add = [[-sc_tar(sv1(1:Nsvnew),:) sc_tar(sv1(1:Nsvnew),:)]; [-sc_nar(sv2(1:Nsvnew),:) sc_nar(sv2(1:Nsvnew),:)]];
end

%-- closest alpha
% score_sv_tar = param.trD.score_SV(param.trD.Y_SV==1,end);
% score_sv_nar = param.trD.score_SV(param.trD.Y_SV==-1,end);
% alpha_tar = param.trD.Alpha(param.trD.Y_SV == 1);
% alpha_nar = param.trD.Alpha(param.trD.Y_SV == -1);
% for ii = 1:Nsv
%     [~,id] = min(abs(score_sv_tar - sc_tar(sv1(ii))));
%     alpha1(ii) = alpha_tar(id);
% end
% for ii = 1:Nsv2
%     [~,id] = min(abs(score_sv_nar - sc_nar(sv2(ii))));
%     alpha2(ii) = alpha_nar(id);
% end

% Alphatemp = [alpha1 alpha2]';
% SV_Y =  [ones(Nsv,1);-ones(Nsv2,1)]; 


Alpha_new_temp = [param.trD.Alpha; Alphatemp];
labels = [param.trD.Y_SV;SV_Y];
SumAlpha = sum(Alpha_new_temp.*labels);
Alpha_new = Alpha_new_temp;

SV_new = [param.trD.SV; SV_add];
Y_SV_new = [param.trD.Y_SV; SV_Y];

while abs(SumAlpha) > 0.0001
    NotMaxAlpha1 = find(Alpha_new <1 & labels == 1);
    NotMaxAlpha2 = find(Alpha_new <1 & labels == -1);
    NNotMaxAlpha = length(NotMaxAlpha1) + length(NotMaxAlpha2);
    deltaAlpha = SumAlpha/NNotMaxAlpha;
    Include1 = Alpha_new(NotMaxAlpha1) >= deltaAlpha & Alpha_new(NotMaxAlpha1) <= 1+deltaAlpha;
    Include2 = Alpha_new(NotMaxAlpha2) <= 1-deltaAlpha & Alpha_new(NotMaxAlpha2) >= -deltaAlpha;
    Alpha_new(NotMaxAlpha1(Include1)) = Alpha_new(NotMaxAlpha1(Include1)) - deltaAlpha;
    Alpha_new(NotMaxAlpha2(Include2)) = Alpha_new(NotMaxAlpha2(Include2)) + deltaAlpha;

%     Alpha_new(Alpha_new < 0) = 0;
%     Alpha_new(Alpha_new > 1) = 1;

    SumAlpha = sum(Alpha_new.*labels);
    %     fprintf('.');
    
    

%     if sum((Alpha_new < 1) &(Alpha_new ~= 0)) == 0
% 
%         NarInd = find(Y_SV_new == -1);
%         TarInd = find(Y_SV_new == 1);
% 
%         Ntar = length(TarInd);
%         Nnar = length(NarInd);
%         NDel = abs(Ntar - Nnar);
% 
%         if Ntar < Nnar
% 
%             NarInd(1:NDel) = [];
%         elseif Ntar > Nar
%             TarInd(1:NDel) = [];
%         end
%         NewInd = sort([TarInd;NarInd]);
% 
%         SumAlpha = sum(Alpha_new(NewInd).*labels(NewInd));
% 
% %         SV_new = SV_new(NewInd,:);
% %         Y_SV_new = Y_SV_new(NewInd,:);
% %         Alpha_new = Alpha_new(NewInd);
% %         score_SV_new = score_SV_new(NewInd);
%     end

end
% fprintf('\n');



zeroalpha = find(Alpha_new == 0);

Alpha_new(zeroalpha) = [];
SV_new(zeroalpha,:) = [];
Y_SV_new(zeroalpha) = [];


score_SV_new = [param.trD.score_SV; sc_add];
score_SV_new(zeroalpha,:) = [];
DelFromNew = find(zeroalpha > size(param.trD.score_SV,1));


% [~,ida] = min(param.trD.Alpha);
Bias_new = mean(Y_SV_new(Alpha_new < 1,:) - sum(mysigmoid(SV_new(Alpha_new > 0,:),SV_new(Alpha_new < 1,:)).*Y_SV_new(Alpha_new > 0,:).*Alpha_new(Alpha_new > 0,:))');

param.trD.Nadded = size(sc_add,1) - length(DelFromNew);
param.trD.Alpha = Alpha_new;
param.trD.Biad =  Bias_new;
param.trD.SV = SV_new;
param.trD.Y_SV = Y_SV_new;
param.trD.score_SV = score_SV_new;
param.trD.score_new = {sc_tar,sc_nar};

end