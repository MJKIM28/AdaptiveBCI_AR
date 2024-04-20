%%
clear all

LAMBDA = 0.02;
%%
fpath = 'E:\1_연구\분석코드\2024_adaptive_BCI\Data\Pilot';
codepath = 'E:\1_연구\분석코드\2024_adaptive_BCI\AdaptiveBCI_AR\0_online_code';


SubNameList = {'Subtest02','Subtest03','Subtest04','Subtest05','Subtest06','Subtest07'};

Nsub = length(SubNameList);

Ntr_tr = 8;
Ntr_con = 15;
Nsess = 6;
Ntr_te = Ntr_con*Nsess;

Nch = 31;
Fs = 500;
Times = -0.2:1/Fs:0.6-1/Fs;
chanloc = readlocs('D:\[1]EEGBCI\[2]Research\Code&Algorithm\BP_channelLocs.locs');
%%
param.DSP.type = 'nf';
param.trD.mode = 'testing';
param.repeat = 10;
param.Stims = 1:4;
param.NumStims = length(param.Stims);
param.Fs = Fs;
param.Baseline = 0.2*param.Fs;
param.Epocline = 0.6*param.Fs;
param.Totalepoc = param.Baseline + param.Epocline;
%%
cd(codepath)
for s = setdiff(1:Nsub,5)
    SubName = SubNameList{s};
    load([fpath,'\Simulate_',num2str(LAMBDA),'param\',SubName])

    load([fpath,'\Test\Epoch\',SubName])

    load([fpath,'\Test\Block\',SubName])
    badch = Block.badch;
    chlist = setdiff(1:Nch,badch);
    chanloc_in = chanloc;
    chanloc_in(badch) = [];
    %%
    upcheck = [UPmodel.upcheck];
    prediction = [UPmodel.prediction];
    target = [UPmodel.target];
    Correctness = prediction == target;
    mdl_init = UPmodel(1).mdl_init;
    DSP_init = UPmodel(1).DSP_init;

     J_ad = []; J_fix = [];
    for t = 1:length(Epoch) - Ntr_con
    
        if t == 1
            DSPweights = DSP_init;
            MDL = mdl_init;
        else
             indices = find(upcheck(1:t-1)==1);
             if ~isempty(indices)
             ind_use = indices(end);           
             DSPweights = UPmodel(ind_use).DSP;
             MDL = UPmodel(ind_use).mdl.mdl;
             else
                DSPweights = DSP_init;
                MDL = mdl_init;
             end

        end
        param.DSP = DSPweights;
        param.trD.mdl = MDL;
        ep = Epoch_condition(Epoch{Ntr_con+t},param);

        [feat] = FeatureExt_DSP(ep,param);
        feat_re  = reshape(feat,param.repeat,param.NumStims,DSPweights.nf,[]);

        %% fisher score (Fixed W, adpative W)
       
        J_ad(:,t) = estimateJ(Epoch{Ntr_con+t},DSPweights.W,DSPweights.theta)';
        J_fix(:,t) = estimateJ(Epoch{Ntr_con+t},DSP_init.W,DSP_init.theta);

    end
    J_allAll{s} = J_ad;
    J_fixAll{s} = J_fix;
    UpCheckAll{s} = upcheck;
    
figure(1);
ii = 1;
ids_temp = find(upcheck);
ids = ids_temp(fix(linspace(1,length(ids_temp),7)));
for i = ids
subplot(Nsub,7,(s-1)*7+ii);
topoplot(UPmodel(i).DSP.W(:,1),chanloc_in); caxis([-1 1]);
ii = ii + 1;
end

figure(2);
    subplot(3,3,s)
    plot(J_ad(1:DSPweights.nf,:)')
    hold on;
    stem(upcheck)
    stem(Correctness*0.5)
   
end




%%
function J= estimateJ(ep_true,DSPweight,theta)

X{1} = ep_true.tar;
X{2} = ep_true.nar;
Nclass = length(X);
[T,C,~] = size(X{1});

S = zeros(C,C,Nclass);
mu = zeros(T,C,Nclass);
for i = 1:Nclass
    N(i) = size(X{i},3);
    
    mu(:,:,i) = mean(X{i},3);
    
    for t = 1:N(i)
        Wd = X{i}(:,:,t) - mu(:,:,i);
        S_i = Wd'*Wd;
        S(:,:,i) = S(:,:,i) + S_i;
    end
    S(:,:,i) = S(:,:,i)/N(i);
end

Sw = sum(S,3);
Sw_reg = (1-theta)*Sw + theta*eye(size(Sw)); % regularize Sw with regularization level theta

mu_all = mean(mu.*permute(repmat(N/sum(N),T,1,C),[1,3,2]),3);

Sb = 0;
for i = 1:Nclass
    Bd = mu(:,:,i) - mu_all;
    Sb = Sb + N(i)*(Bd'*Bd);
end

for comp  = 1:size(DSPweight,2)
w = DSPweight(:,comp);
J(comp) = (w'*Sb*w)/(w'*Sw_reg*w);
end
end