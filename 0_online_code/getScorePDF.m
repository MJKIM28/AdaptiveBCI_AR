function    param = getScorePDF(Feature,label,param)


idt = find(label == 1);
idn = find(label == -1);
idtrain = fix(length(idt)/2);
lbtrain = [label(idt(1:idtrain)); label(idn(1:idtrain*3))];
feattrain = Feature([idt(1:idtrain); idn(1:idtrain*3)],:);

param_temp = param;
param_temp.trD.mode = 'training';
[~,param_temp] = Classification(feattrain,lbtrain,param_temp);

lbtest = label([idt(idtrain+1:end);idn(idtrain*3+1:end)]);
feattest = Feature([idt(idtrain+1:end);idn(idtrain*3+1:end)],:);
[C,sc] = predict(param_temp.trD.mdl,feattest);

sc_train_tar= sc(lbtest == 1,2);
sc_train_nar = sc(lbtest == -1,2);

[px_tar,x_tar] = ksdensity(sc_train_tar,-6:0.1:6);
[px_nar,x_nar] = ksdensity(sc_train_nar,-6:0.1:6);

pdf_score.x = x_tar;
pdf_score.px_tar = px_tar;
pdf_score.px_nar = px_nar;

param.trD.scorePDF = pdf_score;
end