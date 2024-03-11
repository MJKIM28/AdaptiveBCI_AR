function threshold = getthreshold(EP,Feature,label,param)
%%
% idt = find(label == 1);
% idn = find(label == -1);
% idtrain = fix(length(idt)/2);
% lbtrain = [label(idt(1:idtrain)); label(idn(1:idtrain*3))];
% feattrain = Feature([idt(1:idtrain); idn(1:idtrain*3)],:);
% 
% param_temp = param;
% param_temp.trD.mode = 'training';
% [~,param_temp] = ClassificationL(feattrain,lbtrain,param_temp);
% 
% lbtest = label([idt(idtrain+1:end);idn(idtrain*3+1:end)]);
% feattest = Feature([idt(idtrain+1:end);idn(idtrain*3+1:end)],:);
% 
% feat_tar = feattest(lbtest == 1,:);
% feat_nar = feattest(lbtest == -1,:);
% Ntr_temp = ceil(size(feat_tar,1)/param.repeat);
% Nstim_ntar_temp = ceil(size(feat_nar,1)/(param.repeat*Ntr_temp));
% 
% feat_tar_re = reshape(feat_tar',size(feat_tar,2),param.repeat,Ntr_temp);
% feat_nar_re = reshape(feat_nar',size(feat_tar,2),param.repeat,Nstim_ntar_temp,Ntr_temp);
% ZZ = [];
% for tt = 1:Ntr_temp
% Feat = [squeeze(mean(feat_tar_re(:,:,tt),2))'; squeeze(mean(feat_nar_re(:,:,:,tt),2))'];
% [C,sc] = predict(param_temp.trD.mdl,Feat);
% score = sc(:,end);
% ZZ(tt) = max((score - mean(score))/std(score));
% end

% sc_train_tar= sc(lbtest == 1,2);
% sc_train_nar = sc(lbtest == -1,2);
% 
% sc_tar = reshape(sc_train_tar,param.repeat,Ntr_temp);
% sc_tar_mean = mean(reshape(sc_train_tar,param.repeat,Ntr_temp));
% sc_nar = reshape(sc_train_nar,param.repeat*Nstim_ntar_temp,Ntr_temp);
% Dist = [];
% for ii = 1:Ntr_temp
%     Dist(ii) = mahal(sc_tar_mean(:,ii),sc_nar(:,ii));
% end
% threshold = median(Dist);

% threshold = median(ZZ);

%%
global gamma
param_temp = param;
Nnar = param.NumStims - 1;
param_temp.trD.mode = 'training';

EP = Epoch_condition(EP,param_temp);
Ntr = size(EP.tar,3);
Nrp_tr = Ntr/param_temp.repeat_init;
epnar = reshape(EP.nar,param_temp.Totalepoc,param_temp.NumCh,Nnar,param_temp.repeat_init,Nrp_tr);
eptar = reshape(EP.tar,param_temp.Totalepoc,param_temp.NumCh,1,param_temp.repeat_init,Nrp_tr);


s = RandStream('mlfg6331_64');
% figure;
nn = 1;
for tt = 1:4

param_temp.trD.mode = 'training';
param_temp.DSP = rmfield(param_temp.DSP,'W');

% teids = (tt-1)*10+1:tt*10;

EP_train.tar = reshape(eptar(:,:,:,:,setdiff(1:Nrp_tr,tt)),param_temp.Totalepoc,param_temp.NumCh,param_temp.repeat_init*(Nrp_tr-1));
EP_train.nar = reshape(epnar(:,:,:,:,setdiff(1:Nrp_tr,tt)),param_temp.Totalepoc,param_temp.NumCh,param_temp.repeat_init*(Nrp_tr-1)*Nnar);


[feattrain,labeltrain,param_temp] = FeatureExt_DSP(EP_train,param_temp);

gamma = 1/(size(feattrain,2)*var(feattrain(:)));
param_temp.trD.mdl = fitcsvm(feattrain, labeltrain,'KernelFunction','mysigmoid');
param_temp.trD.gamma = gamma;
param_temp.trD.mode      = 'testing';


for ti = 1:100
teids = datasample(s,1:param_temp.repeat_init,param.repeat,'Replace',false);
% teids = 1:10;
EP_test.tar = reshape(eptar(:,:,:,teids,tt),param_temp.Totalepoc,param_temp.NumCh,1,param_temp.repeat);
EP_test.nar = reshape(epnar(:,:,:,teids,tt),param_temp.Totalepoc,param_temp.NumCh,Nnar,param_temp.repeat);
EP_test.nar = cat(3,EP_test.tar,EP_test.nar);

[feattest,labeltest,param_temp] = FeatureExt_DSP(EP_test,param_temp);

gamma = param.trD.gamma;
[C,sc] = predict(param_temp.trD.mdl, feattest);
Nt = size(feattest,1)/(param_temp.repeat*param_temp.NumStims);
score = reshape(sc(:,end),param_temp.repeat,param_temp.NumStims,Nt);

%  score = (score -mean(score(:)))./std(score(:));
[Ntrial,Nstim] = size(param_temp.trD.score);
Dist = [];
for ii = 1:Nstim
    scmean = mean(score(:,ii));
    scother = score(:,setdiff(1:Nstim,ii));
    scothermean = mean(scother(:));
    sgn = sign(scmean - scothermean);
    Dist(ii) = sgn*mahal(scmean,reshape(scother,Ntrial*(Nstim-1),1));

end

[maxval,maxid] = max(Dist);
thretemp(nn) = maxval;
thretempids(nn) = maxid;

% [~,out]= max(mean(param_temp.trD.score))
% clf;
%  gscatter([ones(10,1);2*ones(10,1);3*ones(10,1);4*ones(10,1)],score(:),[ones(10,1);2*ones(10,1);3*ones(10,1);4*ones(10,1)])
% hold on; plot(mean(param_temp.trD.score))
% pause;
nn = nn +1;
end
end
threshold = prctile(thretemp(thretempids == 1),[97.5]);

% 
% feattar = Feature(label==1,:);
% featnar = Feature(label==-1,:);
% Ntr = size(feattar,1);
% Nb = Ntr/param.repeat;
% Nfeat = size(feattar,2);
% feattar_inb = reshape(feattar,param.repeat,Nb,Nfeat);
% featnar_inb = reshape(featnar,Nnar,param.repeat,Nb,Nfeat);
% for ii = 1:Nb
%     traintar = reshape(feattar_inb(:,setdiff(1:Nb,ii),:),param.repeat*(Nb-1),Nfeat);
%     trainnar = reshape(featnar_inb(:,:,setdiff(1:Nb,ii),:),param.repeat*Nnar*(Nb-1),Nfeat);
%     
%     feattrain = [traintar;trainnar];
%     lbtrain = [ones(size(traintar,1),1);-ones(size(trainnar,1),1)];
% 
%     
% 
% 
%     testtar = squeeze(feattar_inb(:,ii,:));
%     testnar = reshape(featnar_inb(:,:,ii,:),Nnar*param.repeat,Nfeat);
% 
%     feattest = [testtar;testnar];
%     [~,sc] = predict(param_temp.trD.mdl,feattest)
% 
% 
% 
% end


end